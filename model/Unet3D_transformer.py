import torch.nn.functional as F
import sys
import os
import torch
import torch.nn as nn
from torch.nn import Linear
from torch_geometric.nn import global_add_pool

# append repo path to sys for convenient imports
sys.path.append(os.path.dirname(os.path.dirname(os.path.realpath(__file__))))

class CryoLigate(nn.Module):
    def __init__(self):
        super().__init__()

        output_shape = (1, 48, 48, 48)
        embedding_size = 768
        
        self.MLP1DTo3D = MLP1DTo3D(embedding_size, output_shape)
        self.UNet3DTransformer = UNet3DTransformer(in_channels=2) 
        self.ReLU = nn.ReLU() 

    def forward(self, data):
        # read data
        low_res_dens, ligand_embedding = data.low_res_dens, data.ligand_embedding
        ligand_embedding = ligand_embedding.squeeze() # remove unnecessary dimensions in the ligands' embeddings
        print("Ligand embedding after squeeze:", ligand_embedding.size())

        # reshape ligand embeddings
        x = self.MLP1DTo3D(ligand_embedding)
        print("After MLP1DTo3D:", x.size())

        # concatenate with low-resolution maps
        x = torch.concat((x, low_res_dens), dim=1)
        print("After concat with low res density:", x.size())

        # feed concatenated data into Unet
        x = self.UNet3DTransformer(x)
        print("After UNet3D:", x.size())
        
        # apply final ReLU
        x = self.ReLU(x)
        print("After final ReLU:", x.size())

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
    
# 3D Convolutional Block
class Conv3DBlock(nn.Module):
    def __init__(self, in_channels, out_channels, bottleneck=False):
        super().__init__()
        self.conv1 = nn.Conv3d(in_channels, out_channels // 2, kernel_size=3, padding=1)
        self.bn1 = nn.BatchNorm3d(out_channels // 2)
        self.conv2 = nn.Conv3d(out_channels // 2, out_channels, kernel_size=3, padding=1)
        self.bn2 = nn.BatchNorm3d(out_channels)
        self.relu = nn.ReLU()
        self.bottleneck = bottleneck
        if not bottleneck:
            self.pooling = nn.MaxPool3d(kernel_size=2, stride=2)
    
    def forward(self, x):
        x = self.relu(self.bn1(self.conv1(x)))
        x = self.relu(self.bn2(self.conv2(x)))
        return (self.pooling(x) if not self.bottleneck else x), x

# Transformer Encoder Block with Attention
class Transformer3D(nn.Module):
    def __init__(self, embed_dim, num_heads, mlp_dim):
        super().__init__()
        self.norm1 = nn.LayerNorm(embed_dim)
        self.attn = nn.MultiheadAttention(embed_dim, num_heads)
        self.norm2 = nn.LayerNorm(embed_dim)
        self.mlp = nn.Sequential(
            nn.Linear(embed_dim, mlp_dim),
            nn.ReLU(),
            nn.Linear(mlp_dim, embed_dim)
        )
    
    def forward(self, x):
        b, c, d, h, w = x.shape
        x = x.view(b, c, -1).permute(2, 0, 1)  # Flatten spatial dims, shape: (D*H*W, B, C)
        x = self.norm1(x + self.attn(x, x, x)[0])
        x = self.norm2(x + self.mlp(x))
        x = x.permute(1, 2, 0).view(b, c, d, h, w)
        return x

# Squeeze-and-Excitation Attention
class SEAttention(nn.Module):
    def __init__(self, channels, reduction=16):
        super().__init__()
        self.fc1 = nn.Linear(channels, channels // reduction)
        self.relu = nn.ReLU()
        self.fc2 = nn.Linear(channels // reduction, channels)
        self.sigmoid = nn.Sigmoid()
    
    def forward(self, x):
        b, c, d, h, w = x.shape
        y = x.view(b, c, -1).mean(dim=-1)  # Global Average Pooling
        y = self.relu(self.fc1(y))
        y = self.sigmoid(self.fc2(y)).view(b, c, 1, 1, 1)
        return x * y

# UpSampling Block with Attention
class UpConv3DBlock(nn.Module):
    def __init__(self, in_channels, res_channels=0, last_layer=False):
        super().__init__()
        self.upconv = nn.ConvTranspose3d(in_channels, in_channels, kernel_size=2, stride=2)
        self.conv1 = nn.Conv3d(in_channels + res_channels, in_channels // 2, kernel_size=3, padding=1)
        self.conv2 = nn.Conv3d(in_channels // 2, in_channels // 2, kernel_size=3, padding=1)
        self.bn = nn.BatchNorm3d(in_channels // 2)
        self.relu = nn.ReLU()
        self.attention = SEAttention(in_channels // 2)
        self.last_layer = last_layer
        if last_layer:
            self.conv3 = nn.Conv3d(in_channels // 2, 1, kernel_size=1)  # Single-channel output
    
    def forward(self, x, res=None):
        x = self.upconv(x)
        if res is not None:
            x = torch.cat((x, res), dim=1)
        x = self.relu(self.bn(self.conv1(x)))
        x = self.relu(self.bn(self.conv2(x)))
        x = self.attention(x)
        return self.conv3(x) if self.last_layer else x

# 3D UNet with Transformer and Attention
class UNet3DTransformer(nn.Module):
    def __init__(self, in_channels=1, level_channels=[64, 128, 256], bottleneck_channel=512):
        super().__init__()
        self.encoder1 = Conv3DBlock(in_channels, level_channels[0])
        self.encoder2 = Conv3DBlock(level_channels[0], level_channels[1])
        self.encoder3 = Conv3DBlock(level_channels[1], level_channels[2])
        self.bottleneck = Conv3DBlock(level_channels[2], bottleneck_channel, bottleneck=True)
        self.transformer = Transformer3D(bottleneck_channel, num_heads=8, mlp_dim=1024)
        self.decoder3 = UpConv3DBlock(bottleneck_channel, level_channels[2])
        self.decoder2 = UpConv3DBlock(level_channels[2], level_channels[1])
        self.decoder1 = UpConv3DBlock(level_channels[1], level_channels[0], last_layer=True)
    
    def forward(self, x):
        residual = x  # Save input for residual learning
        x, res1 = self.encoder1(x)
        x, res2 = self.encoder2(x)
        x, res3 = self.encoder3(x)
        x, _ = self.bottleneck(x)
        x = self.transformer(x)  # Transformer bottleneck
        x = self.decoder3(x, res3)
        x = self.decoder2(x, res2)
        x = self.decoder1(x, res1)
        return x + residual  # Residual connection for sharpening

# Frequency-based Fourier Loss
# def fft_loss(pred, target):
#     pred_fft = torch.fft.fftn(pred, dim=(-3, -2, -1))
#     target_fft = torch.fft.fftn(target, dim=(-3, -2, -1))
#     return F.l1_loss(torch.abs(pred_fft), torch.abs(target_fft))

# Example Loss Usage
# criterion = lambda pred, target: 0.7 * CustomLoss()(pred, target) + 0.3 * fft_loss(pred, target)
