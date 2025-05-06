import torch
import sys
import pickle
import mrcfile
import joblib
import os
import numpy as np
from torch_geometric.data import Batch, Data

# append repo path to sys for convenient imports
sys.path.append(os.path.dirname(os.path.dirname(os.path.realpath(__file__))))

from data_generation.generate_dataset import normalize_density_minmax
from utils import delete_extension_from_filename, extract_filename_from_full_path


def generate_pyg(ligand_embedding_path, low_res_density_path, save_path):
    """
    Generates a .pyg file for testing.

    Args:
        ligand_embedding_path - full path to the file with ligand embedding
        low_res_density_path - full path to the file with input low resolution density 
        save_path - full path to the output .pyg file 
    """
    
    # read ligand embedding data
    ligand_embedding = torch.load(ligand_embedding_path)

    # read and normalize input low resolution density
    low_res_density = mrcfile.read(low_res_density_path)
    low_res_density_normalized = normalize_density_minmax(low_res_density)
    low_res_dens = torch.tensor(low_res_density_normalized).unsqueeze(0).unsqueeze(0)

    # save data to a .pyg file
    data = Data(
        ligand_embedding=ligand_embedding,
        low_res_dens=low_res_dens
    )   
    torch.save(data, save_path)



if __name__ == "__main__":
    
    # path to the folder with the experimnetal maps/ligands/embeddings etc
    # maps processed by the network are also stored in this folder
    experimental_data = os.path.join(os.getcwd(), "experimental_data")
    
    # list with the downsampled densities names
    low_res_densities_names = [
        "6ply_Ligand_Map_filtered_gaus_0.7_newscale_newbox.mrc",
        "8hnd_Ligand_Map_filtered_gaus_0.7_newscale_newbox.mrc",
        "8w8s_Ligand_Map_filtered_gaus_0.7_newscale_newbox.mrc",
        "8txz_Ligand_Map_filtered_gaus_0.7_newscale_newbox.mrc",
        "6w6e_Ligand_Map_filtered_gaus_0.7_newscale_newbox.mrc",
        "8hfl_Ligand_Map_filtered_gaus_0.7_newscale_newbox.mrc",
        "8xv2_Ligand_Map_filtered_gaus_0.7_newscale_newbox.mrc"
    ]

    # list with the corresponding embeddings names
    # NOTE: generate embeddings and fill this list before run the code!
    embeddings_names = [
        "",
        "",
        "",
        "",
        "",
        "",
        ""
    ]


    # low_res_densities_names = [
    #     "6ply_Ligand_Map_filtered_gaus_0.7_newscale_newbox.mrc",
    # ]

    # # list with the corresponding embeddings names
    # # NOTE: generate embeddings and fill this list before run the code!
    # embeddings_names = [
    #     "",
    # ]


    # specify the model and corresponding epoch (on which the model was saved) to process data
    epoch = 79
    model_path = os.path.join(
        os.getcwd(),
        "model_save_folder",
        "SCUnet3D_without_Ligand_embeddings_20.0L1_0.1SSIM_loss_Norm_minmax_maps_Forward_model_bad_nconfs3_to_Good_res2.0_Batchsize_32_lr_5.0e-04_wd_1.0e-04_20250502_102045",
        "model",
        f"model_{epoch}.pkl"
    )  
    model_name = f"SCUnet3D_20.0L1_0.1SSIM_Batchsize_32_Epoch{epoch}" # shorter name for the model to use in the processed map names 
    
    # load the model
    device = "cuda" if torch.cuda.is_available() else "cpu"
    model = joblib.load(model_path)
    model.to(device)
    model.eval() 


    # generate maps for all complexes from the lists above
    for i in range(len(low_res_densities_names)):
        # construct full paths to the necessary files and generate a .pyg file for the network
        low_res_density_name = low_res_densities_names[i]
        embedding_name = embeddings_names[i]
        ligand_embedding_path = os.path.join(experimental_data, embedding_name)
        low_res_density_path = os.path.join(experimental_data, low_res_density_name)
        save_path = delete_extension_from_filename(low_res_density_path) + "_with_embedding.pyg"
        generate_pyg(ligand_embedding_path, low_res_density_path, save_path)

        # read an input pyg file and generate prediction
        data = torch.load(save_path)
        data = data.to(device)
        with torch.no_grad():
            pred = model(data)
            pred = pred[0][0]
        pred = pred.cpu()
        pred = pred.detach().numpy()
        print(f"min: {np.min(pred)}")
        print(f"max: {np.max(pred)}")
        print(f"mean: {np.mean(pred)}")
        print(f"median: {np.median(pred)}")
        
        # extract voxel size and origin from the input low resolution map 
        # we need this information to properly save predicted map
        low_res_dens_origin = np.array([0, 0, 0])
        low_res_voxel_size = (0.5, 0.5, 0.5)
        with mrcfile.open(low_res_density_path, "r") as low_res_dens_file:
            low_res_dens_origin = low_res_dens_file.header.origin
            low_res_voxel_size = low_res_dens_file.voxel_size

        # save predicted map with the proper origin and voxel size
        predicted_map_path = os.path.join(
                        experimental_data,
                        f"{delete_extension_from_filename(extract_filename_from_full_path(low_res_density_path))}_resolved_with_{model_name}.mrc"
                        )        
        with mrcfile.new(predicted_map_path, overwrite=True) as predicted_file:
            predicted_file.set_data(pred)
            predicted_file.header.origin = low_res_dens_origin
            predicted_file.voxel_size = low_res_voxel_size
    
