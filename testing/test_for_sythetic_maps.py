import torch
import sys
import pickle
import mrcfile
import joblib
import os
import numpy as np
from sklearn.model_selection import KFold
import pandas as pd

# append repo path to sys for convenient imports
sys.path.append(os.path.dirname(os.path.dirname(os.path.realpath(__file__))))

from data_generation.generate_dataset import NetworkDataset, PLIDataLoader

device = "cuda" if torch.cuda.is_available() else "cpu"



data_root = '/proj/berzelius-2022-haloi/users/x_nanha'
toy_dir = os.path.join(data_root, 'PDBBind_Zenodo_6408497')
toy_df = pd.read_csv(os.path.join(toy_dir, "Test_for_sythetic_map_train.csv"))
# toy_df = pd.read_csv(os.path.join(toy_dir, "PDB_IDs_with_rdkit_length_less_than_24A_nbonds_less_than_10.csv"))
# n_splits = 10
# kfold = KFold(n_splits=n_splits, shuffle=True, random_state=42)
# kfold_split = list(kfold.split(toy_df))
# k = 0
# train_idx, val_idx = kfold_split[k]
# train_df = toy_df.iloc[train_idx]
valid_df = toy_df

# print(train_df)

num_process = 20
create_dataset = False
base_ligand_filename = "_ligand.pdb"
base_ligand_embedding_filename = "_ligand_embedding.pyg"
base_label_filename = "_ligand_res_2.0_gridpsace_0.5_nbox_48_size_24A_label.mrc"
base_low_res_density_filename = "_nconfs3_genmode_gnina_docking_boxextens1.0_res4.0_delprob0.0_low_resolution_forward_model.mrc"
clarifying_data_name = "_forward_model_nconfs3_to_good_res2.0_norm_minmax_with_embeddings"
is_dataset_log = True
dataset_log_path = os.path.join(os.getcwd(), "network_dataset_main_logs")

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
valid_loader = PLIDataLoader(valid_set, batch_size=1, shuffle=False, num_workers=4)

data = next(iter(valid_loader))
data = data.to(device)
# print(data)

# pdb_id = "1i43"
# example = f"/proj/berzelius-2022-haloi/users/x_nanha/PDBBind_Zenodo_6408497/{pdb_id}/{pdb_id}_ligand_forward_model_nconfs3_to_good_res2.0_norm_minmax_with_embeddings,pyg"


epoch = 79
model_path = os.path.join(
    os.getcwd(),
    "model_save_folder",
    "SCUnet3D_without_Ligand_embeddings_20.0L1_0.1SSIM_loss_Norm_minmax_maps_Forward_model_bad_nconfs3_to_Good_res2.0_Batchsize_32_lr_5.0e-04_wd_1.0e-04_20250502_102045",
    "model",
    f"model_{epoch}.pkl"
)  
model_name = f"SCUnet3D_20.0L1_0.1SSIM_Batchsize_32_Epoch{epoch}" # shorter name for the model to use in the processed map names 

# # load the model
model = joblib.load(model_path)
model.to(device)
model.eval()

with torch.no_grad():
    pred = model(data)
    pred = pred[0][0]
pred = pred.cpu()
pred = pred.detach().numpy()
mrcfile.write(f'testing/1t31_generated.mrc', pred, overwrite=True)