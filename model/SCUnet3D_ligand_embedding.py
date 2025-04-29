# Copyright (C) 2024  Hong Cao, Jiahua He, Tao Li, Sheng-You Huang and Huazhong University of Science
# Copyright (C) 2022  Kai Zhang (cskaizhang@gmail.com, https://cszn.github.io/).

# All rights reserved.

# 3D version of Swin-Conv-UNet (https://arxiv.org/pdf/2203.13278.pdf)

# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#     http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.

import sys
import os
import torch
import numpy as np
import torch.nn as nn
from einops import rearrange 
from einops.layers.torch import Rearrange
from timm.models.layers import trunc_normal_, DropPath

sys.path.append(os.path.dirname(os.path.dirname(os.path.realpath(__file__))))

from model.frn import FilterResponseNorm3d


class WMSA(nn.Module):
    def __init__(self, input_dim, output_dim, head_dim, window_size, type):
        super(WMSA, self).__init__()
        self.input_dim = input_dim
        self.output_dim = output_dim
        self.head_dim = head_dim 
        self.scale = self.head_dim ** -0.5
        self.n_heads = input_dim//head_dim
        self.window_size = window_size
        self.type=type
        self.embedding_layer = nn.Linear(self.input_dim, 3*self.input_dim, bias=True) 

        self.relative_position_params = nn.Parameter(torch.zeros((2 * window_size - 1)*(2 * window_size -1)*(2 * window_size -1), self.n_heads))

        self.linear = nn.Linear(self.input_dim, self.output_dim)

        trunc_normal_(self.relative_position_params, std=.02)
        self.relative_position_params = torch.nn.Parameter(self.relative_position_params.view(2*window_size-1, 2*window_size-1, 2*window_size-1, self.n_heads).transpose(2,3).transpose(1,2).transpose(0,1))

        # Modified by Jiahua He 
        cord = torch.tensor(np.array([[i, j, k] for i in range(self.window_size) for j in range(self.window_size) for k in range(self.window_size)]))
        self.relation = cord[:, None, :] - cord[None, :, :] + self.window_size -1

    def generate_mask(self, h, w, d, p, shift):
        attn_mask = torch.zeros(h, w, d, p, p, p, p, p, p, dtype=torch.bool, device=self.relative_position_params.device)
        if self.type == 'W':
            return attn_mask

        s = p - shift
        attn_mask[-1, :,  :, :s, :,  :,  s:, :,  :] = True
        attn_mask[-1, :,  :, s:, :,  :,  :s, :,  :] = True
        attn_mask[:, -1,  :, :, :s,  :,  :, s:,  :] = True
        attn_mask[:, -1,  :, :, s:,  :,  :, :s,  :] = True
        attn_mask[:,  :, -1, :,  :, :s,  :,  :, s:] = True
        attn_mask[:,  :, -1, :,  :, s:,  :,  :, :s] = True
        attn_mask = rearrange(attn_mask, 'w1 w2 w3 p1 p2 p3 p4 p5 p6 -> 1 1 (w1 w2 w3) (p1 p2 p3) (p4 p5 p6)')

        return attn_mask

    def forward(self, x):
        if self.type!='W': x = torch.roll(x, shifts=(-(self.window_size//2), -(self.window_size//2), -(self.window_size//2)), dims=(1,2,3))
        x = rearrange(x, 'b (w1 p1) (w2 p2) (w3 p3) c -> b w1 w2 w3 p1 p2 p3 c', p1=self.window_size, p2=self.window_size, p3=self.window_size)
        h_windows = x.size(1)
        w_windows = x.size(2)
        d_windows = x.size(3)
        assert h_windows == w_windows == d_windows

        x = rearrange(x, 'b w1 w2 w3 p1 p2 p3 c -> b (w1 w2 w3) (p1 p2 p3) c', p1=self.window_size, p2=self.window_size, p3=self.window_size)

        qkv = self.embedding_layer(x)

        q, k, v = rearrange(qkv, 'b nw np (threeh c) -> threeh b nw np c', c=self.head_dim).chunk(3, dim=0)
        sim = torch.einsum('hbwpc,hbwqc->hbwpq', q, k) * self.scale

        sim = sim + rearrange(self.relative_position_params[:, self.relation[:,:,0].long(), self.relation[:,:,1].long(), self.relation[:,:,2].long()], 'h p q -> h 1 1 p q')

        if self.type != 'W':
            attn_mask = self.generate_mask(h_windows, w_windows, d_windows, self.window_size, shift=self.window_size//2)
            sim = sim.masked_fill_(attn_mask, float("-inf"))

        probs = nn.functional.softmax(sim, dim=-1)
        output = torch.einsum('hbwij,hbwjc->hbwic', probs, v)
        output = rearrange(output, 'h b w p c -> b w p (h c)')
        output = self.linear(output)
        output = rearrange(output, 'b (w1 w2 w3) (p1 p2 p3) c -> b (w1 p1) (w2 p2) (w3 p3) c', w1=h_windows, w2=w_windows, w3=d_windows, p1=self.window_size, p2=self.window_size, p3=self.window_size)

        if self.type!='W': output = torch.roll(output, shifts=(self.window_size//2, self.window_size//2, self.window_size//2), dims=(1,2,3))
        return output

class Block(nn.Module):
    def __init__(self, input_dim, output_dim, head_dim, window_size, drop_path, type='W', input_resolution=None):
        super(Block, self).__init__()
        self.input_dim = input_dim
        self.output_dim = output_dim
        assert type in ['W', 'SW']
        self.type = type
        if input_resolution <= window_size:
            self.type = 'W'

        print("Block Initial Type: {}, drop_path_rate:{:.6f}".format(self.type, drop_path))
        self.ln1 = nn.LayerNorm(input_dim)
        self.msa = WMSA(input_dim, input_dim, head_dim, window_size, self.type)
        self.drop_path = DropPath(drop_path) if drop_path > 0. else nn.Identity()
        self.ln2 = nn.LayerNorm(input_dim)
        self.mlp = nn.Sequential(
            nn.Linear(input_dim, 4 * input_dim),
            nn.GELU(),
            nn.Linear(4 * input_dim, output_dim),
        )

    def forward(self, x):
        x = self.ln1(x)
        print("After LayerNorm-7asdsa:", x.size())

        x = self.msa(x)

        x = self.drop_path(x)

        # x = x + self.drop_path(self.msa(self.ln1(x)))

        x = x + self.drop_path(self.mlp(self.ln2(x)))
        return x

class ConvTransBlock(nn.Module):
    def __init__(self, conv_dim, trans_dim, head_dim, window_size, drop_path, type='W', input_resolution=None):
        super(ConvTransBlock, self).__init__()
        self.conv_dim = conv_dim
        self.trans_dim = trans_dim
        self.head_dim = head_dim
        self.window_size = window_size
        self.drop_path = drop_path
        self.type = type
        self.input_resolution = input_resolution

        assert self.type in ['W', 'SW']
        if self.input_resolution <= self.window_size:
            self.type = 'W'

        self.trans_block = Block(self.trans_dim, self.trans_dim, self.head_dim, self.window_size, self.drop_path, self.type, self.input_resolution)
        # This is Conv3d-2 as per in the paper Suppl data 10
        self.conv1_1 = nn.Conv3d(self.conv_dim+self.trans_dim, self.conv_dim+self.trans_dim, 1, 1, 0, bias=True)
        self.conv1_2 = nn.Conv3d(self.conv_dim+self.trans_dim, self.conv_dim+self.trans_dim, 1, 1, 0, bias=True)

        self.conv_block = nn.Sequential(
                nn.Conv3d(self.conv_dim, self.conv_dim, 3, 1, 1, bias=False),
                FilterResponseNorm3d(self.conv_dim),
                nn.Conv3d(self.conv_dim, self.conv_dim, 3, 1, 1, bias=False),
                FilterResponseNorm3d(self.conv_dim),
                )

    def forward(self, x):
        conv1_1 = self.conv1_1(x)
        print("After ConvTransBlock Conv3d-2:", conv1_1.size())
        
        conv_x, trans_x = torch.split(conv1_1, (self.conv_dim, self.trans_dim), dim=1)
        print("After ConvTransBlock conv_x trans_x:", conv_x.size(), trans_x.size())
        
        conv_x = self.conv_block(conv_x) + conv_x
        print("After Going through the RConv block:", conv_x.size())

        trans_x = Rearrange('b c h w d -> b h w d c')(trans_x)
        print("After Rearranging trans_x:", trans_x.size())

        trans_x = self.trans_block(trans_x)
        trans_x = Rearrange('b h w d c -> b c h w d')(trans_x)
        res = self.conv1_2(torch.cat((conv_x, trans_x), dim=1))
        x = x + res

        return x


class MLP1DTo3D(nn.Module):
    def __init__(self, input_size, output_shape):
        super(MLP1DTo3D, self).__init__()
        self.output_shape = output_shape  # Shape you want to reshape to (A, B, C)
        
        level_channels=[256, 256*16, 256*32] 
        hidden_size_1, hidden_size_2, hidden_size_3 = level_channels[0], level_channels[1], level_channels[2]
        
        # Define the MLP layers
        self.fc1 = nn.Linear(input_size, hidden_size_1)
        self.fc2 = nn.Linear(hidden_size_1, hidden_size_2)
        self.fc3 = nn.Linear(hidden_size_2, hidden_size_3)
        self.fc4 = nn.Linear(hidden_size_3, output_shape[0] * output_shape[1] * output_shape[2] * output_shape[3])
    
    def forward(self, x):
        x = torch.relu(self.fc1(x))
        x = torch.relu(self.fc2(x))
        x = torch.relu(self.fc3(x))
        
        # Output the desired number of elements to reshape into 3D
        x = self.fc4(x)
        
        # Reshape the output to 3D
        x = x.view(-1, *self.output_shape)  # Reshapes output to (batch_size, A, B, C)
        
        return x
    
class SCUNet(nn.Module):

    def __init__(self, in_nc=2, config=[2,2,2,2,2,2,2], dim=32, drop_path_rate=0.2, input_resolution=48, head_dim=16, window_size=3, n_classes=1):
        super(SCUNet, self).__init__()
        self.config = config
        self.dim = dim
        self.head_dim = head_dim
        self.window_size = window_size

        output_shape = (1, 48, 48, 48)
        embedding_size = 768
        
        self.MLP1DTo3D = MLP1DTo3D(embedding_size, output_shape)

        dpr = [x.item() for x in torch.linspace(0, drop_path_rate, sum(config))]

        self.m_head = [nn.Conv3d(in_nc, dim, 3, 1, 1, bias=False)]

        begin = 0
        self.m_down1 = [ConvTransBlock(dim//2, dim//2, self.head_dim, self.window_size, dpr[i+begin], 'W' if not i%2 else 'SW', input_resolution) 
                      for i in range(config[0])] + \
                      [nn.Conv3d(dim, 2*dim, 2, 2, 0, bias=False)]

        begin += config[0]
        self.m_down2 = [ConvTransBlock(dim, dim, self.head_dim, self.window_size, dpr[i+begin], 'W' if not i%2 else 'SW', input_resolution//2)
                      for i in range(config[1])] + \
                      [nn.Conv3d(2*dim, 4*dim, 2, 2, 0, bias=False)]

        begin += config[1]
        self.m_down3 = [ConvTransBlock(2*dim, 2*dim, self.head_dim, self.window_size, dpr[i+begin], 'W' if not i%2 else 'SW',input_resolution//4)
                      for i in range(config[2])] + \
                      [nn.Conv3d(4*dim, 8*dim, 2, 2, 0, bias=False)]

        begin += config[2]
        self.m_body = [ConvTransBlock(4*dim, 4*dim, self.head_dim, self.window_size, dpr[i+begin], 'W' if not i%2 else 'SW', input_resolution//8)
                    for i in range(config[3])]

        begin += config[3]
        self.m_up3 = [nn.ConvTranspose3d(8*dim, 4*dim, 2, 2, 0, bias=False),] + \
                      [ConvTransBlock(2*dim, 2*dim, self.head_dim, self.window_size, dpr[i+begin], 'W' if not i%2 else 'SW',input_resolution//4)
                      for i in range(config[4])]
                      
        begin += config[4]
        self.m_up2 = [nn.ConvTranspose3d(4*dim, 2*dim, 2, 2, 0, bias=False),] + \
                      [ConvTransBlock(dim, dim, self.head_dim, self.window_size, dpr[i+begin], 'W' if not i%2 else 'SW', input_resolution//2)
                      for i in range(config[5])]
                      
        begin += config[5]
        self.m_up1 = [nn.ConvTranspose3d(2*dim, dim, 2, 2, 0, bias=False),] + \
                    [ConvTransBlock(dim//2, dim//2, self.head_dim, self.window_size, dpr[i+begin], 'W' if not i%2 else 'SW', input_resolution) 
                      for i in range(config[6])]

        self.m_tail = [nn.Conv3d(dim, n_classes, 3, 1, 1, bias=False)]

        self.m_head = nn.Sequential(*self.m_head)
        self.m_down1 = nn.Sequential(*self.m_down1)
        self.m_down2 = nn.Sequential(*self.m_down2)
        self.m_down3 = nn.Sequential(*self.m_down3)
        self.m_body = nn.Sequential(*self.m_body)
        self.m_up3 = nn.Sequential(*self.m_up3)
        self.m_up2 = nn.Sequential(*self.m_up2)
        self.m_up1 = nn.Sequential(*self.m_up1)
        self.m_tail = nn.Sequential(*self.m_tail)  

    def forward(self, data):
        low_res_dens, ligand_embedding = data.low_res_dens, data.ligand_embedding
        print("Ligand embedding before squeeze shape:", ligand_embedding.size())
        print("low_res_dens shape:", low_res_dens.size())

        ligand_embedding = ligand_embedding.squeeze() # remove unnecessary dimensions in the ligands' embeddings
        print("Ligand embedding shape:", ligand_embedding.size())

        # reshape ligand embeddings
        ligand_embedding = self.MLP1DTo3D(ligand_embedding)
        print("Ligand embedding after MLP1DTo3D shape:", ligand_embedding.size())

        # concatenate with low-resolution maps
        x0 = torch.concat((ligand_embedding, low_res_dens), dim=1)
        print("Ligand embedding concat with Low res density shape:", x0.size())

        x1 = self.m_head(x0)
        print("After Conv3d-1:", x1.size())
        x2 = self.m_down1(x1)
        x3 = self.m_down2(x2)
        x4 = self.m_down3(x3)
        x = self.m_body(x4)
        x = self.m_up3(x+x4)
        x = self.m_up2(x+x3)
        x = self.m_up1(x+x2)
        x = self.m_tail(x+x1)

        return x

    def _init_weights(self, m):
        if isinstance(m, nn.Linear):
            trunc_normal_(m.weight, std=.02)
            if m.bias is not None:
                nn.init.constant_(m.bias, 0)
        elif isinstance(m, nn.LayerNorm):
            nn.init.constant_(m.bias, 0)
            nn.init.constant_(m.weight, 1.0)


# if __name__ == '__main__':

#     # torch.cuda.empty_cache()
#     net = SCUNet()
#     print("Num params: ", sum(p.numel() for p in net.parameters()))
#     x = torch.randn((1, 1, 48, 48, 48))
#     x = net(x)
#     # print(x.shape)