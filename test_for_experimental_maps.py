import torch
import pickle
import mrcfile
import joblib
import os
import numpy as np
from generate_dataset import generate_network_data
from torch_geometric.data import Batch, Data
from utils import delete_extension_from_filename, extract_filename_from_full_path


def generate_pyg(low_res_density_path, save_path):
    
    low_res_density = mrcfile.read(low_res_density_path)

    low_res_density = normalize_density_minmax(low_res_density)
    # low_res_density = normalize_density_with_percentile(low_res_density, percentile=0.99)

    low_res_dens = torch.tensor(low_res_density).unsqueeze(0).unsqueeze(0)

    
    data = Data(
        low_res_dens=low_res_dens
    )

    torch.save(data, save_path)



if __name__ == "__main__":
    experimental_data = os.path.join(os.getcwd(), "experimental_data")

    low_res_densities_names = [
        "6ply_Ligand_Map_filtered_gaus_0.7_newscale_newbox.mrc",
        "8hnd_Ligand_Map_filtered_gaus_0.7_newscale_newbox.mrc",
        "8w8s_Ligand_Map_filtered_gaus_0.7_newscale_newbox.mrc",
        "8txz_Ligand_Map_filtered_gaus_0.7_newscale_newbox.mrc",
        "6w6e_Ligand_Map_filtered_gaus_0.7_newscale_newbox.mrc",
        "8hfl_Ligand_Map_filtered_gaus_0.7_newscale_newbox.mrc",
        "8xv2_Ligand_Map_filtered_gaus_0.7_newscale_newbox.mrc"
    ]
    for low_res_density_name in low_res_densities_names:
        low_res_density_path = os.path.join(experimental_data, low_res_density_name)
        save_path = delete_extension_from_filename(low_res_density_path) + ".pyg"
        generate_pyg(low_res_density_path, save_path)


        epoch = 253

        model_name = f"L2_with_ReLU_k7_epoch{epoch}"
        model_path = os.path.join(
            os.getcwd(),
            "model",
            "k_7_only_final_unet_with_final_ReLU_L2loss_norm_minmax_maps_forward_model_bad_nconfs3_to_good_res2.0_batchsize_64_hidsize_256_levels_256_lr_5e-4_wd_1e-5_20250106_161358",
            "model",
            f"model_{epoch}.pkl"
            )  


        device = "cuda" if torch.cuda.is_available() else "cpu"

        data = torch.load(save_path)
        data = data.to(device)

        model = joblib.load(model_path)
        model.to(device)
        model.eval() 

        with torch.no_grad():
            pred = model(data)
            pred = pred[0][0]

        pred = pred.cpu()
        pred = pred.detach().numpy()

        print(f"min: {np.min(pred)}")
        print(f"max: {np.max(pred)}")
        print(f"mean: {np.mean(pred)}")
        print(f"median: {np.median(pred)}")

        low_res_dens_origin = np.array([0, 0, 0])
        low_res_voxel_size = (0.5, 0.5, 0.5)

        with mrcfile.open(low_res_density_path, "r") as low_res_dens_file:
            low_res_dens_origin = low_res_dens_file.header.origin
            low_res_voxel_size = low_res_dens_file.voxel_size

        generated_map_path = os.path.join(
                        experimental_data,
                        f"{delete_extension_from_filename(extract_filename_from_full_path(low_res_density_path))}_resolved_with_{model_name}.mrc"
                        )
        
        with mrcfile.new(generated_map_path, overwrite=True) as generated_file:
            generated_file.set_data(pred)
            generated_file.header.origin = low_res_dens_origin
            generated_file.voxel_size = low_res_voxel_size
    
