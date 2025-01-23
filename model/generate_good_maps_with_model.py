import os
import sys
import joblib
import pandas as pd
os.environ['CUDA_VISIBLE_DEVICES'] = "0"
import torch
torch.cuda.empty_cache()
import numpy as np
from sklearn.model_selection import KFold
import mrcfile

# append repo path to sys for convenient imports
sys.path.append(os.path.dirname(os.path.dirname(os.path.realpath(__file__))))

from data_generation.generate_dataset import NetworkDataset, PLIDataLoader

from utils import *

if __name__ == "__main__":

    # path to the database
    data_root = '/proj/berzelius-2022-haloi/users/x_elima'
    toy_dir = os.path.join(data_root, 'PDBBind_Zenodo_6408497')

     # specify the model and corresponding epoch (on which the model was saved) to process data
    epoch = 28
    model_name = f"k0_unet_with_embeddings_combined_loss_{epoch}"
    model_path = os.path.join(
            os.getcwd(),
            "model",
            "k_0_unet_with_embeddings_combined_loss_norm_minmax_maps_forward_model_bad_nconfs3_to_good_res2.0_batchsize_32_lr_5e-4_wd_1e-4_20250117_132748",
            "model",
            f"model_{epoch}.pkl"
            )  

    # repeat cross-valiation scenario to properly split dataset into train and validation
    toy_df = pd.read_csv(os.path.join(toy_dir, "PDB_IDs_with_rdkit_length_less_than_24A_nbonds_less_than_10.csv"))
    n_splits = 10
    kfold = KFold(n_splits=n_splits, shuffle=True, random_state=42)
    kfold_split = list(kfold.split(toy_df))
    k = 0 # NOTE: specify the index of the split corresponding to the model
    train_idx, val_idx = kfold_split[k]
    train_df = toy_df.iloc[train_idx]
    valid_df = toy_df.iloc[val_idx]

    # extract one example from the training set and one from the validation to generate maps
    # NOTE: you can specify whichever indices you want, not necessarily the same
    train_complex_id = 333
    val_complex_id = 333
    train_df = train_df.iloc[train_complex_id: train_complex_id + 1]
    valid_df = valid_df.iloc[val_complex_id: val_complex_id + 1]
    complex_name_train = train_df.iloc[0]["pdbid"]
    complex_name_val = valid_df.iloc[0]["pdbid"]
    print(f"Generate map with {model_name} for complex {complex_name_train} from train")
    print(f"Generate map with {model_name} for complex {complex_name_val} from validation")

    # generate input network data for train and validation complexes
    num_process = 20
    create_dataset = False
    base_ligand_filename = "_ligand.pdb"
    base_ligand_embedding_filename = "_ligand_embedding.pyg"
    base_label_filename = "_ligand_res_2.0_gridpsace_0.5_nbox_48_size_24A_label.mrc"
    base_low_res_density_filename = "_nconfs3_genmode_gnina_docking_boxextens1.0_res4.0_delprob0.0_low_resolution_forward_model.mrc"
    clarifying_data_name = "_forward_model_nconfs3_to_good_res2.0_norm_minmax_with_embeddings"
    is_dataset_log = True
    dataset_log_path = os.path.join(os.getcwd(), "network_dataset_main_logs")
    train_set = NetworkDataset(
                toy_dir, 
                train_df, 
                num_process=num_process,
                create=create_dataset,
                base_ligand_filename=base_ligand_filename,
                base_ligand_embedding_filename=base_ligand_embedding_filename,
                base_label_filename=base_label_filename,
                base_low_res_density_filename=base_low_res_density_filename,
                clarifying_data_name=clarifying_data_name,
                is_log=is_dataset_log,
                log_path=dataset_log_path,
            )
    valid_set = NetworkDataset(
                toy_dir, 
                valid_df, 
                num_process=num_process,
                create=create_dataset,
                base_ligand_filename=base_ligand_filename,
                base_ligand_embedding_filename=base_ligand_embedding_filename,
                base_label_filename=base_label_filename,
                base_low_res_density_filename=base_low_res_density_filename,
                clarifying_data_name=clarifying_data_name,
                is_log=is_dataset_log,
                log_path=dataset_log_path,
            )
    train_loader = PLIDataLoader(train_set, batch_size=1, shuffle=False, drop_last=True)
    valid_loader = PLIDataLoader(valid_set, batch_size=1, shuffle=False, num_workers=4)

    # load the model
    device = "cuda" if torch.cuda.is_available() else "cpu"
    print("Device: ", device)
    model = joblib.load(model_path)
    model.to(device)
    model.eval()

    # generate map for the complex from training set
    for data in train_loader:
        data = data.to(device)
        print("Train Data size", data.size())
        with torch.no_grad():
            pred = model(data)
            pred = pred[0][0]
        pred = pred.cpu()
        pred = pred.detach().numpy()
        print("Train PRED SHAPE", pred.shape)

        # extract voxel size and origin from the input low resolution map 
        # we need this information to properly save predicted map
        low_res_dens_path = os.path.join(
                        toy_dir,
                        complex_name_train,
                        complex_name_train + base_low_res_density_filename
                        )        
        low_res_dens_origin = np.array([0, 0, 0])
        low_res_voxel_size = (0.5, 0.5, 0.5)
        with mrcfile.open(low_res_dens_path, "r") as low_res_dens_file:
            low_res_dens_origin = low_res_dens_file.header.origin
            low_res_voxel_size = low_res_dens_file.voxel_size

        # save generatred map with the proper origin and voxel size
        generated_map_path = os.path.join(
                        os.getcwd(),
                        "maps_generated_with_model",
                        complex_name_train + "_train_" + model_name + "_generated_map.mrc"
                        )        
        with mrcfile.new(generated_map_path, overwrite=True) as generated_file:
            generated_file.set_data(pred)
            generated_file.header.origin = low_res_dens_origin
            generated_file.voxel_size = low_res_voxel_size


    # repeat for the complex from training set
    for data in valid_loader:
        data = data.to(device)
        print("Val Data size", data.size())
        with torch.no_grad():
            pred = model(data)
            pred = pred[0][0]

        pred = pred.cpu()
        pred = pred.detach().numpy()

        print("VAL PRED SHAPE", pred.shape)

        low_res_dens_path =  os.path.join(
                        toy_dir,
                        complex_name_val,
                        complex_name_val + base_low_res_density_filename
                        )
        
        low_res_dens_origin = np.array([0, 0, 0])
        low_res_voxel_size = (0.5, 0.5, 0.5)
        

        with mrcfile.open(low_res_dens_path, "r") as low_res_dens_file:
            low_res_dens_origin = low_res_dens_file.header.origin
            low_res_voxel_size = low_res_dens_file.voxel_size

        generated_map_path = os.path.join(
                        os.getcwd(),
                        "maps_generated_with_model",
                        complex_name_val + "_val_" + model_name + "_generated_map.mrc"
                        )
        
        with mrcfile.new(generated_map_path, overwrite=True) as generated_file:
            generated_file.set_data(pred)
            generated_file.header.origin = low_res_dens_origin
            generated_file.voxel_size = low_res_voxel_size
