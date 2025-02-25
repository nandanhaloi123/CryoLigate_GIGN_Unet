import torch
import torch.nn as nn
import torch.nn.functional as F

class CryoLigate(nn.Module):
    def __init__(self):
        super().__init__()

        output_shape = (1, 48, 48, 48)
        embedding_size = 768

        self.MLP1DTo3D = MLP1DTo3D(embedding_size, output_shape)
        self.UNet3DTransformer = UNet3DTransformer(in_channels=2) 
        self.ReLU = nn.ReLU() 

    def forward(self, data):
        low_res_dens, ligand_embedding = data.low_res_dens, data.ligand_embedding
        ligand_embedding = ligand_embedding.squeeze()  # Remove unnecessary dimensions

        x = self.MLP1DTo3D(ligand_embedding)
        x = torch.cat((x, low_res_dens), dim=1)  # Concatenate along channel dimension

        x = self.UNet3DTransformer(x)
        x = self.ReLU(x)

        return x


class MLP1DTo3D(nn.Module):
    def __init__(self, input_size, output_shape):
        super().__init__()
        self.output_shape = output_shape  

        level_channels = [256, 256 * 16, 256 * 32] 
        hidden_size_1, hidden_size_2, hidden_size_3 = level_channels[0], level_channels[1], level_channels[2]

        self.fc1 = nn.Linear(input_size, hidden_size_1)
        self.fc2 = nn.Linear(hidden_size_1, hidden_size_2)
        self.fc3 = nn.Linear(hidden_size_2, hidden_size_3)
        self.fc4 = nn.Linear(hidden_size_3, output_shape[0] * output_shape[1] * output_shape[2] * output_shape[3])

    def forward(self, x):
        x = torch.relu(self.fc1(x))
        x = torch.relu(self.fc2(x))
        x = torch.relu(self.fc3(x))
        x = self.fc4(x)
        x = x.view(-1, *self.output_shape)  # Reshape into 3D
        return x


class Conv3DBlock(nn.Module):
    """3D Convolutional Block with Transformer."""
    def __init__(self, in_channels, out_channels, apply_transformer=True, bottleneck=False):
        super().__init__()
        self.conv1 = nn.Conv3d(in_channels, out_channels // 2, kernel_size=3, padding=1)
        self.bn1 = nn.BatchNorm3d(out_channels // 2)
        self.conv2 = nn.Conv3d(out_channels // 2, out_channels, kernel_size=3, padding=1)
        self.bn2 = nn.BatchNorm3d(out_channels)
        self.relu = nn.ReLU()
        self.bottleneck = bottleneck
        self.apply_transformer = apply_transformer

        if not bottleneck:
            self.pooling = nn.MaxPool3d(kernel_size=2, stride=2)

        # Transformer for feature refinement
        if apply_transformer:
            self.transformer = Transformer3D(out_channels, num_heads=8, mlp_dim=1024)

    def forward(self, x):
        x = self.relu(self.bn1(self.conv1(x)))
        x = self.relu(self.bn2(self.conv2(x)))

        if self.apply_transformer:
            x = self.transformer(x)

        return (self.pooling(x) if not self.bottleneck else x), x


class Transformer3D(nn.Module):
    """Transformer block for 3D feature refinement."""
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
        x = x.view(b, c, -1).permute(2, 0, 1)  # Flatten spatial dims -> (D*H*W, B, C)
        x = self.norm1(x + self.attn(x, x, x)[0])
        x = self.norm2(x + self.mlp(x))
        x = x.permute(1, 2, 0).view(b, c, d, h, w)
        return x


class UpConv3DBlock(nn.Module):
    """3D UpSampling Block with Transformer."""
    def __init__(self, in_channels, res_channels=0, apply_transformer=True, last_layer=False):
        super().__init__()
        self.upconv = nn.ConvTranspose3d(in_channels, in_channels, kernel_size=2, stride=2)
        self.conv1 = nn.Conv3d(in_channels + res_channels, in_channels // 2, kernel_size=3, padding=1)
        self.conv2 = nn.Conv3d(in_channels // 2, in_channels // 2, kernel_size=3, padding=1)
        self.bn = nn.BatchNorm3d(in_channels // 2)
        self.relu = nn.ReLU()
        self.last_layer = last_layer

        # Transformer block for decoder layers
        if apply_transformer:
            self.transformer = Transformer3D(in_channels // 2, num_heads=8, mlp_dim=1024)

        if last_layer:
            self.conv3 = nn.Conv3d(in_channels // 2, 1, kernel_size=1)  # Single-channel output

    def forward(self, x, res=None):
        x = self.upconv(x)
        if res is not None:
            x = torch.cat((x, res), dim=1)
        x = self.relu(self.bn(self.conv1(x)))
        x = self.relu(self.bn(self.conv2(x)))
        
        if self.apply_transformer:
            x = self.transformer(x)

        return self.conv3(x) if self.last_layer else x


class UNet3DTransformer(nn.Module):
    """3D UNet with Transformer at each encoder and decoder stage."""
    def __init__(self, in_channels=1, level_channels=[64, 128, 256], bottleneck_channel=512):
        super().__init__()
        self.encoder1 = Conv3DBlock(in_channels, level_channels[0], apply_transformer=True)
        self.encoder2 = Conv3DBlock(level_channels[0], level_channels[1], apply_transformer=True)
        self.encoder3 = Conv3DBlock(level_channels[1], level_channels[2], apply_transformer=True)
        self.bottleneck = Conv3DBlock(level_channels[2], bottleneck_channel, apply_transformer=True, bottleneck=True)
        
        # Decoders
        self.decoder3 = UpConv3DBlock(bottleneck_channel, level_channels[2], apply_transformer=False)
        self.decoder2 = UpConv3DBlock(level_channels[2], level_channels[1], apply_transformer=False)
        self.decoder1 = UpConv3DBlock(level_channels[1], level_channels[0], apply_transformer=False, last_layer=True)

    def forward(self, x):
        residual = x  # Save input for residual connection

        x, res1 = self.encoder1(x)
        x, res2 = self.encoder2(x)
        x, res3 = self.encoder3(x)
        x, _ = self.bottleneck(x)

        x = self.decoder3(x, res3)
        x = self.decoder2(x, res2)
        x = self.decoder1(x, res1)

        return x + residual  # Residual connection for sharpening
