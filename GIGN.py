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

        output_shape = (1, 32, 32, 32)
        
        self.MLP1DTo3D = MLP1DTo3D(hidden_dim, output_shape)
        self.UNet3D_low_res_dens = UNet3D(in_channels=1, num_classes=1) # NOTE TODO: add/remove additional Unet for low resolution maps
        self.UNet3D = UNet3D(in_channels=1, num_classes=1) #NOTE TODO: change in_channels=2 if use low resolution maps
        self.ReLU = nn.ReLU() # NOTE TODO: add/remove final ReLU
        #self.fc = FC(hidden_dim, hidden_dim, 3, 0.1, 1)

    def forward(self, data):
        #NOTE TODO: return low_res_dens
        x, edge_index_intra, edge_index_inter, pos, low_res_dens = \
        data.x, data.edge_index_intra, data.edge_index_inter, data.pos, data.low_res_dens
        # x, edge_index_intra, edge_index_inter, pos = \
        # data.x, data.edge_index_intra, data.edge_index_inter, data.pos

        # # NOTE TODO: return Graph NN
        # x = self.lin_node(x)
        # x = self.gconv1(x, edge_index_intra, edge_index_inter, pos)
        # x = self.gconv2(x, edge_index_intra, edge_index_inter, pos)
        # x = self.gconv3(x, edge_index_intra, edge_index_inter, pos)
        # # print("Before global_add_pool:", x.size())
        # x = global_add_pool(x, data.batch)
        # # print("After global_add_pool:", x.size())
        # x = self.MLP1DTo3D(x)
        # # print("After MLP1DTo3D:", x.size())

        # # NOTE TODO: add/remove additional Unet for low resolution maps
        # low_res_dens_after_unet = self.UNet3D_low_res_dens(low_res_dens)
        # print("Low res dens after Unet:", low_res_dens_after_unet.size())

        # # NOTE TODO: return concatenation
        # x = torch.concat((x, low_res_dens_after_unet), dim=1)
        # print("After concat with low res density:", x.size())

        # NOTE TODO: return Unet3D
        x = self.UNet3D(low_res_dens)
        # print("After UNet3D:", x.size())
        
        # # NOTE TODO: add/remove final RelU
        # x = self.ReLU(x)
        # print("After final ReLU:", x.size())

        # x = self.fc(x)
        # print("After FC:", x.size())
        return x
        #return x.view(-1)
    
    
class MLP1DTo3D(nn.Module):
    def __init__(self, input_size, output_shape):
        super(MLP1DTo3D, self).__init__()
        self.output_shape = output_shape  # Shape you want to reshape to (A, B, C)
        
        level_channels=[256, 256*16, 256*32] 
        # level_channels=[256*2, 256*16, 256*32] # NOTE TODO: change levels here
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

# Example usage:
# input_size = 64  # 1D input size
# hidden_size = 128  # Hidden layer size
# output_shape = (4, 4, 4)  # Desired 3D shape (A, B, C)




# class FC(nn.Module):
#     def __init__(self, d_graph_layer, d_FC_layer, n_FC_layer, dropout, n_tasks):
#         super(FC, self).__init__()
#         self.d_graph_layer = d_graph_layer
#         self.d_FC_layer = d_FC_layer
#         self.n_FC_layer = n_FC_layer
#         self.dropout = dropout
#         self.predict = nn.ModuleList()
#         for j in range(self.n_FC_layer):
#             if j == 0:
#                 self.predict.append(nn.Linear(self.d_graph_layer, self.d_FC_layer))
#                 self.predict.append(nn.Dropout(self.dropout))
#                 self.predict.append(nn.LeakyReLU())
#                 self.predict.append(nn.BatchNorm1d(d_FC_layer))
#             if j == self.n_FC_layer - 1:
#                 self.predict.append(nn.Linear(self.d_FC_layer, n_tasks))
#             else:
#                 self.predict.append(nn.Linear(self.d_FC_layer, self.d_FC_layer))
#                 self.predict.append(nn.Dropout(self.dropout))
#                 self.predict.append(nn.LeakyReLU())
#                 self.predict.append(nn.BatchNorm1d(d_FC_layer))

#     def forward(self, h):
#         for layer in self.predict:
#             h = layer(h)

#         return h
    
#device = torch.device('cuda:0')
#model = GIGN(35, 256).to(device)
#print("Num params: ", sum(p.numel() for p in model.parameters()))
# print(model)
