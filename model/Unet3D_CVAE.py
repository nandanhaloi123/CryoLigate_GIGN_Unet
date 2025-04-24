import torch.nn.functional as F
import sys
import os
import torch
import torch.nn as nn
from torch.nn import Linear
from torch_geometric.nn import global_add_pool

# append repo path to sys for convenient imports
sys.path.append(os.path.dirname(os.path.dirname(os.path.realpath(__file__))))

class CryoLigateCVAE(nn.Module):
    def __init__(self):
        super().__init__()
        
        Ligand_output_shape_after_MLP = (1, 48, 48, 48)
        embedding_size = 768
        latent_dim = 256  # Latent space dimension
        
        self.MLP1DTo3D = MLP1DTo3D(embedding_size, Ligand_output_shape_after_MLP)
        self.UNet3DTransformerCVAE = UNet3DTransformerCVAE(in_channels=2, latent_dim=latent_dim, embedding_dim=embedding_size)
        self.ReLU = nn.ReLU()
    
    def forward(self, data):
        low_res_dens, ligand_embedding = data.low_res_dens, data.ligand_embedding
        print("Ligand embedding before squeeze shape:", ligand_embedding.size())
        ligand_embedding = ligand_embedding.squeeze()
        print("Ligand embedding shape:", ligand_embedding.size())

        x = self.MLP1DTo3D(ligand_embedding)
        print("Ligand embedding after MLP1DTo3D shape:", x.size())
        x = torch.concat((x, low_res_dens), dim=1)
        print("Ligand embedding concat with Low res density shape:", x.size())
        
        x, mu, logvar, actual_noise, predicted_noise = self.UNet3DTransformerCVAE(x, ligand_embedding)
        x = self.ReLU(x)
        return x, mu, logvar, actual_noise, predicted_noise


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


class Conv3DBlock(nn.Module):
    """
    The basic block for double 3x3x3 convolutions in the analysis path
    -- __init__()
    :param in_channels -> number of input channels
    :param out_channels -> desired number of output channels
    :param bottleneck -> specifies the bottlneck block
    -- forward()
    :param input -> input Tensor to be convolved
    :return -> Tensor
    """

    def __init__(self, in_channels, out_channels, bottleneck = False) -> None:
        super(Conv3DBlock, self).__init__()
        self.conv1 = nn.Conv3d(in_channels= in_channels, out_channels=out_channels//2, kernel_size=(3,3,3), padding=1)
        self.bn1 = nn.BatchNorm3d(num_features=out_channels//2)
        self.conv2 = nn.Conv3d(in_channels= out_channels//2, out_channels=out_channels, kernel_size=(3,3,3), padding=1)
        self.bn2 = nn.BatchNorm3d(num_features=out_channels)
        self.relu = nn.ReLU()
        self.bottleneck = bottleneck
        if not bottleneck:
            self.pooling = nn.MaxPool3d(kernel_size=(2,2,2), stride=2)

    
    def forward(self, input):
        print("Encoder block 1 shape:", input.size())
        res = self.relu(self.bn1(self.conv1(input)))
        print("Encoder block 2 shape:", res.size())
        res = self.relu(self.bn2(self.conv2(res)))
        print("Encoder block 3 shape:", res.size())
        out = None
        if not self.bottleneck:
            out = self.pooling(res)
            print("Encoder block after pooling shape:", res.size())
        else:
            out = res
        return out, res

class UpConv3DBlock(nn.Module):
    """
    The basic block for upsampling followed by double 3x3x3 convolutions in the synthesis path
    -- __init__()
    :param in_channels -> number of input channels
    :param out_channels -> number of residual connections' channels to be concatenated
    :param last_layer -> specifies the last output layer
    :param num_classes -> specifies the number of output channels for dispirate classes
    -- forward()
    :param input -> input Tensor
    :param residual -> residual connection to be concatenated with input
    :return -> Tensor
    """

    def __init__(self, in_channels, res_channels=0, last_layer=False, num_classes=None) -> None:
        super(UpConv3DBlock, self).__init__()
        assert (last_layer==False and num_classes==None) or (last_layer==True and num_classes!=None), 'Invalid arguments'
        self.upconv1 = nn.ConvTranspose3d(in_channels=in_channels, out_channels=in_channels, kernel_size=(2, 2, 2), stride=2)
        self.relu = nn.ReLU()
        self.bn = nn.BatchNorm3d(num_features=in_channels//2)
        self.conv1 = nn.Conv3d(in_channels=in_channels+res_channels, out_channels=in_channels//2, kernel_size=(3,3,3), padding=(1,1,1))
        self.conv2 = nn.Conv3d(in_channels=in_channels//2, out_channels=in_channels//2, kernel_size=(3,3,3), padding=(1,1,1))
        self.last_layer = last_layer
        if last_layer:
            self.conv3 = nn.Conv3d(in_channels=in_channels//2, out_channels=num_classes, kernel_size=(1,1,1))
            
        
    def forward(self, input, residual=None):
        print("Decoder block before upconv shape:", input.size())
        out = self.upconv1(input)
        print("Decoder block 1 shape:", out.size())
        if residual!=None: out = torch.cat((out, residual), 1)
        print("Decoder block 1 --after concat shape:", out.size())
        out = self.relu(self.bn(self.conv1(out)))
        print("Decoder block 2 shape:", out.size())
        out = self.relu(self.bn(self.conv2(out)))
        print("Decoder block 3 shape:", out.size())
        if self.last_layer: out = self.conv3(out)
        return out

# 3D UNet with Transformer and Attention
class UNet3DTransformerCVAE(nn.Module):
    def __init__(self, in_channels=2, level_channels=[64, 128, 256], bottleneck_channel=512, latent_dim=256, embedding_dim=768):
        super().__init__()
        self.encoder1 = Conv3DBlock(in_channels, level_channels[0])
        self.encoder2 = Conv3DBlock(level_channels[0], level_channels[1])
        self.encoder3 = Conv3DBlock(level_channels[1], level_channels[2])
        self.bottleneck = Conv3DBlock(level_channels[2], bottleneck_channel, bottleneck=True)
        
        # Conditioning mechanism for ligand embedding
        self.embed_fc = nn.Linear(embedding_dim, latent_dim)
        self.fc_mu = nn.Linear(bottleneck_channel * 6 * 6 * 6 + latent_dim, latent_dim)
        self.fc_logvar = nn.Linear(bottleneck_channel * 6 * 6 * 6 + latent_dim, latent_dim)
        self.fc_decode = nn.Linear(latent_dim + latent_dim, bottleneck_channel * 6 * 6 * 6)
        
        self.transformer = Transformer3D(bottleneck_channel, num_heads=8, mlp_dim=1024)
        
        self.decoder3 = UpConv3DBlock(bottleneck_channel, level_channels[2])
        self.decoder2 = UpConv3DBlock(level_channels[2], level_channels[1])
        self.decoder1 = UpConv3DBlock(level_channels[1], level_channels[0], last_layer=True)
    
    def reparameterize(self, mu, logvar):
        std = torch.exp(0.5 * logvar)
        eps = torch.randn_like(std)
        return mu + eps * std
    
    def forward(self, x, ligand_embedding):
        x, res1 = self.encoder1(x)
        x, res2 = self.encoder2(x)
        x, res3 = self.encoder3(x)
        x, _ = self.bottleneck(x)
        
        x_flat = x.view(x.size(0), -1)
        ligand_embed = self.embed_fc(ligand_embedding)
        
        x_cat = torch.cat([x_flat, ligand_embed], dim=-1)
        mu = self.fc_mu(x_cat)
        logvar = self.fc_logvar(x_cat)
        z = self.reparameterize(mu, logvar)
        
        z_cat = torch.cat([z, ligand_embed], dim=-1)
        x = self.fc_decode(z_cat).view(x.size(0), -1, 6, 6, 6)
        x = self.transformer(x)
        
        x = self.decoder3(x, res3)
        x = self.decoder2(x, res2)
        x = self.decoder1(x, res1)
        return x, mu, logvar