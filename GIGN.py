# %%
import torch
import torch.nn as nn
from torch.nn import Linear
from torch_geometric.nn import global_add_pool
from HIL import HIL
from Unet3D import *

class GIGN(nn.Module):
    def __init__(self, node_dim, hidden_dim):
        super().__init__()
        self.lin_node = nn.Sequential(Linear(node_dim, hidden_dim), nn.SiLU())
    
        self.gconv1 = HIL(hidden_dim, hidden_dim)
        self.gconv2 = HIL(hidden_dim, hidden_dim)
        self.gconv3 = HIL(hidden_dim, hidden_dim)

        output_shape = (1, 48, 48, 48)
        embedding_size = 768
        
        self.MLP1DTo3D = MLP1DTo3D(embedding_size, output_shape)
        self.UNet3D = UNet3D(in_channels=2, num_classes=1) 
        self.ReLU = nn.ReLU() 
        #self.fc = FC(hidden_dim, hidden_dim, 3, 0.1, 1)

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
        x = self.UNet3D(x)
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
