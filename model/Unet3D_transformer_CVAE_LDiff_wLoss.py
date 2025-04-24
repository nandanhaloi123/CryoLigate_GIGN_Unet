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


class LatentDiffusion(nn.Module):
    def __init__(self, latent_dim, num_steps=1000, beta_start=0.0001, beta_end=0.02):
        super().__init__()
        self.num_steps = num_steps
        self.betas = torch.linspace(beta_start, beta_end, num_steps)
        self.alphas = 1 - self.betas
        self.alphas_cumprod = torch.cumprod(self.alphas, dim=0)

        self.noise_predictor = nn.Sequential(
            nn.Linear(latent_dim, 512),
            nn.ReLU(),
            nn.Linear(512, latent_dim)
        )

    def forward(self, z, t=None):
        batch_size = z.shape[0]
        device = z.device  # Get the device of input tensor `z`
        
        # Ensure alphas_cumprod is on the same device as the input z
        self.alphas_cumprod = self.alphas_cumprod.to(device)

        if t is None:
            t = torch.randint(0, self.num_steps, (batch_size,), device=device)  # Ensure `t` is on the same device as `z`

        noise = torch.randn_like(z, device=device)  # Ensure `noise` is on the same device as `z`

        # Apply noise
        alpha_t = self.alphas_cumprod[t].view(-1, 1)  # `alpha_t` is already on the correct device
        z_noisy = torch.sqrt(alpha_t) * z + torch.sqrt(1 - alpha_t) * noise

        # Predict noise
        predicted_noise = self.noise_predictor(z_noisy)

        # Remove noise
        z_denoised = (z_noisy - predicted_noise) / torch.sqrt(alpha_t)

        return z_denoised, noise, predicted_noise
    
# 3D UNet with Transformer and Attention
class UNet3DTransformerCVAE(nn.Module):
    def __init__(self, in_channels=2, level_channels=[64, 128, 256], bottleneck_channel=512, latent_dim=256, embedding_dim=768):
        super().__init__()
        self.encoder1 = Conv3DBlock(in_channels, level_channels[0])
        self.encoder2 = Conv3DBlock(level_channels[0], level_channels[1])
        self.encoder3 = Conv3DBlock(level_channels[1], level_channels[2])
        self.bottleneck = Conv3DBlock(level_channels[2], bottleneck_channel, bottleneck=True)
        
        # NEED A NETWORK WITH NON LINEARITY HERE?
        Ligand_output_shape_after_MLP = (1, 6, 6, 6)
        self.MLP1DTo3D = MLP1DTo3D(embedding_dim, Ligand_output_shape_after_MLP)

        # Conditioning mechanism for ligand embedding
        # self.embed_fc = nn.Linear(embedding_dim, latent_dim) 

        self.fc_mu = nn.Linear((bottleneck_channel+1) * 6 * 6 * 6, latent_dim)
        self.fc_logvar = nn.Linear((bottleneck_channel+1) * 6 * 6 * 6, latent_dim)
        
        self.latent_diffusion = LatentDiffusion(latent_dim)  # Add Latent Diffusion here

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
        print("Entering bottleneck:", x.size())
        x, _ = self.bottleneck(x)

        lig_mlp = self.MLP1DTo3D(ligand_embedding)
        print("Ligand embedding after MLP1DTo3D shape:", lig_mlp.size())

        x = torch.concat((x, lig_mlp), dim=1)
        print("Ligand embedding concat with last layer of bottleneck:", x.size())
        
        x_flat = x.view(x.size(0), -1)
        print("Flattened shape:", x_flat.size()) # DO WE NEED THE FLATTENING HERE?

        # ligand_embed = self.embed_fc(ligand_embedding)
        # print("Ligand embedding before adding to flattened vector:", ligand_embed.size()) #NEED SOME WAYS OF CONCATENATING HERE.

        # x_cat = torch.cat([x_flat, ligand_embed], dim=-1)
        # print("Flattened + Ligand embedding shape:", x_cat.size())
        mu = self.fc_mu(x_flat)
        print("Mu shape:", mu.size())
        logvar = self.fc_logvar(x_flat)
        print("logvar shape:", logvar.size())
        z = self.reparameterize(mu, logvar)
        print("After sampling Z shape:", z.size())
        
        # Apply Latent Diffusion
        z, actual_noise, predicted_noise  = self.latent_diffusion(z)

        print("After diffusion Z shape:", z.size())
        print("After diffusion actual noise shape:", actual_noise.size())
        print("After diffusion predicted noise shape:", predicted_noise.size())

        # z_cat = torch.cat([z, ligand_embed], dim=-1)
        x = self.fc_decode(z).view(x.size(0), -1, 6, 6, 6)
        x = self.transformer(x)
        
        x = self.decoder3(x, res3)
        x = self.decoder2(x, res2)
        x = self.decoder1(x, res1)
        return x, mu, logvar, actual_noise, predicted_noise

