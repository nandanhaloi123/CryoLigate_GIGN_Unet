import os
import joblib
os.environ['CUDA_VISIBLE_DEVICES'] = "0"
import torch
torch.cuda.empty_cache()
import torch.nn as nn
import torch.optim as optim
import pandas as pd
from utils import AverageMeter
from datetime import datetime, timezone
from GIGN import GIGN
from dataset_GIGN import GraphDataset, PLIDataLoader
from config.config_dict import Config
from log.train_logger import TrainLogger
import numpy as np
from utils import *
from sklearn.metrics import mean_squared_error
from sklearn.model_selection import KFold
import mrcfile



if __name__ == "__main__":
    data_root = '/proj/berzelius-2022-haloi/users/x_elima'

    toy_dir = os.path.join(data_root, 'PDBBind_Zenodo_6408497')


    epoch = 60

    model_name = f"L2_with_ReLU_epoch{epoch}"
    model_path = os.path.join(
            os.getcwd(),
            "model",
            "only_final_unet_with_final_ReLU_L2loss_norm_minmax_maps_forward_model_bad_nconfs3_to_good_res2.0_batchsize_64_hidsize_256_levels_256_lr_5e-4_wd_1e-5_20250103_170201",
            "model",
            f"model_{epoch}.pkl"
            )  

    
    ###################################
    toy_df = pd.read_csv(os.path.join(toy_dir, "PDB_IDs_with_rdkit_length_less_than_24A.csv")).sample(frac=1., random_state=123)
    split_idx = int(0.9 * len(toy_df))
    train_df = toy_df.iloc[:split_idx]
    valid_df = toy_df.iloc[split_idx:]

    train_complex_id = 255
    val_complex_id = 255
    train_df = train_df.iloc[train_complex_id: train_complex_id + 1]
    valid_df = valid_df.iloc[val_complex_id: val_complex_id + 1]


    complex_name_train = train_df.iloc[0]["pdbid"]
    complex_name_val = valid_df.iloc[0]["pdbid"]

    print(f"Generate maps with {model_name} for complex {complex_name_train} from train")
    print(f"Generate maps with {model_name} for complex {complex_name_val} from validation")

    graph_type = "Graph_GIGN"
    dis_threshold = 5
    num_process = 20
    create_dataset = False
    # base_corr_check_filename = "_Nandan_graphs_with_bad_maps.txt"
    # base_corr_check_filename = "_corr0.6_passed_complexes.txt"
    base_corr_check_filename = "_complexes_temporary_bad_maps_only.txt"
    base_label_filename = "_ligand_res_2.0_gridpsace_0.5_nbox_48_size_24A_label.mrc"
    base_low_res_density_filename = "_nconfs3_genmode_gnina_docking_boxextens1.0_res4.0_delprob0.0_low_resolution_forward_model.mrc"
    clarifying_graph_name = "_forward_model_bad_nconfs3_to_good_res2.0_norm_minmax"
    is_dataset_log = True
    dataset_log_path = os.path.join(os.getcwd(), "dataset_GIGN_main_logs")

    train_set = GraphDataset(
                toy_dir, 
                train_df, 
                dis_threshold=dis_threshold,
                graph_type=graph_type, 
                num_process=num_process,
                create=create_dataset,
                base_corr_check_filename=base_corr_check_filename,
                base_label_filename=base_label_filename,
                base_low_res_density_filename=base_low_res_density_filename,
                clarifying_graph_name=clarifying_graph_name,
                is_log=is_dataset_log,
                log_path=dataset_log_path,
            )
    valid_set = GraphDataset(
                toy_dir, 
                valid_df, 
                dis_threshold=dis_threshold,
                graph_type=graph_type, 
                num_process=num_process,
                create=create_dataset,
                base_corr_check_filename=base_corr_check_filename,
                base_label_filename=base_label_filename,
                base_low_res_density_filename=base_low_res_density_filename,
                clarifying_graph_name=clarifying_graph_name,
                is_log=is_dataset_log,
                log_path=dataset_log_path,
            )

    train_loader = PLIDataLoader(train_set, batch_size=1, shuffle=False, drop_last=True)
    valid_loader = PLIDataLoader(valid_set, batch_size=1, shuffle=False, num_workers=4)


    device = "cuda" if torch.cuda.is_available() else "cpu"
    model = joblib.load(model_path)
    model.to(device)
    model.eval()


    for data in train_loader:
        data = data.to(device)
        print("Train Data size", data.size())
        with torch.no_grad():
            pred = model(data)
            pred = pred[0][0]

        pred = pred.cpu()
        pred = pred.detach().numpy()

        print("Train PRED SHAPE", pred.shape)

        low_res_dens_path =  os.path.join(
                        toy_dir,
                        complex_name_train,
                        complex_name_train + base_low_res_density_filename
                        )
        
        low_res_dens_origin = np.array([0, 0, 0])
        low_res_voxel_size = (0.5, 0.5, 0.5)

        with mrcfile.open(low_res_dens_path, "r") as low_res_dens_file:
            low_res_dens_origin = low_res_dens_file.header.origin
            low_res_voxel_size = low_res_dens_file.voxel_size

        generated_map_path = os.path.join(
                        os.getcwd(),
                        "maps_generated_with_model",
                        complex_name_train + "_train_" + model_name + "_generated_map.mrc"
                        )
        
        with mrcfile.new(generated_map_path, overwrite=True) as generated_file:
            generated_file.set_data(pred)
            generated_file.header.origin = low_res_dens_origin
            generated_file.voxel_size = low_res_voxel_size



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

    ###########################################

    #toy_df = pd.read_csv(os.path.join(toy_dir, "PDB_IDs_with_rdkit_length_less_than_16A_succ_gnina.csv"))
    # k_folds = 5
    # kfold = KFold(n_splits=k_folds, shuffle=True, random_state=42)
    # for fold, (train_idx, val_idx) in enumerate(kfold.split(toy_df)):
    #     train_complex_id = train_idx[0]
    #     val_complex_id = val_idx[0]

    #     train_df = toy_df.iloc[train_complex_id: train_complex_id + 1]
    #     valid_df = toy_df.iloc[val_complex_id: val_complex_id + 1]


    #     complex_name_train = train_df.iloc[0]["pdbid"]
    #     complex_name_val = valid_df.iloc[0]["pdbid"]

    #     print(f"Generate maps with {model_name} for complex {complex_name_train} from train")
    #     print(f"Generate maps with {model_name} for complex {complex_name_val} from validation")

    #     graph_type = "Graph_GIGN"
    #     dis_threshold = 5
    #     num_process = 20
    #     create_dataset = False
    #     base_corr_check_filename = "_Nandan_graphs_with_bad_maps.txt"
    #     # base_corr_check_filename = "_Nandan_graphs.txt"
    #     # base_corr_check_filename = "_corr0.6_passed_complexes.txt"
    #     base_label_filename = "_ligand_res1.0_boxed_16A.mrc"
    #     base_low_res_density_filename = "_nconfs10_genmodedocking_res3.5_nbox16_threshcorr0.3_delprob0.2_low_resolution_forward_model.mrc"
    #     is_dataset_log = True
    #     dataset_log_path = os.path.join(os.getcwd(), "dataset_GIGN_main_logs")

    #     train_set = GraphDataset(
    #                 toy_dir, 
    #                 train_df, 
    #                 dis_threshold=dis_threshold,
    #                 graph_type=graph_type, 
    #                 num_process=num_process,
    #                 create=create_dataset,
    #                 base_corr_check_filename=base_corr_check_filename,
    #                 base_label_filename=base_label_filename,
    #                 base_low_res_density_filename=base_low_res_density_filename,
    #                 is_log=is_dataset_log,
    #                 log_path=dataset_log_path,
    #             )
        
    #     valid_set = GraphDataset(
    #                 toy_dir, 
    #                 valid_df, 
    #                 dis_threshold=dis_threshold,
    #                 graph_type=graph_type, 
    #                 num_process=num_process,
    #                 create=create_dataset,
    #                 base_corr_check_filename=base_corr_check_filename,
    #                 base_label_filename=base_label_filename,
    #                 base_low_res_density_filename=base_low_res_density_filename,
    #                 is_log=is_dataset_log,
    #                 log_path=dataset_log_path,
    #             )
        
    #     train_loader = PLIDataLoader(train_set, batch_size=1, shuffle=False, num_workers=4)
    #     valid_loader = PLIDataLoader(valid_set, batch_size=1, shuffle=False, num_workers=4)


    #     device = "cuda" if torch.cuda.is_available() else "cpu"
    #     model = joblib.load(model_path)
    #     model.to(device)
    #     model.eval()


    #     for data in train_loader:
    #         data = data.to(device)
    #         print("Train Data size", data.size())
    #         with torch.no_grad():
    #             pred = model(data)
    #             pred = pred[0][0]

    #         pred = pred.cpu()
    #         pred = pred.detach().numpy()

    #         print("Train PRED SHAPE", pred.shape)

    #         low_res_dens_path =  os.path.join(
    #                         toy_dir,
    #                         complex_name_train,
    #                         complex_name_train + base_low_res_density_filename
    #                         )
            
    #         low_res_dens_origin = np.array([0, 0, 0])

    #         with mrcfile.open(low_res_dens_path, "r") as low_res_dens_file:
    #             low_res_dens_origin = low_res_dens_file.header.origin

    #         generated_map_path = os.path.join(
    #                         toy_dir,
    #                         complex_name_train,
    #                         complex_name_train + "_train_" + model_name + "_generated_map.mrc"
    #                         )
            
    #         with mrcfile.new(generated_map_path, overwrite=True) as generated_file:
    #             generated_file.set_data(pred)
    #             generated_file.header.origin = low_res_dens_origin
    #             generated_file.voxel_size = (0.5, 0.5, 0.5)


    
    #     for data in valid_loader:
    #         data = data.to(device)
    #         print("Val Data size", data.size())
    #         with torch.no_grad():
    #             pred = model(data)
    #             pred = pred[0][0]

    #         pred = pred.cpu()
    #         pred = pred.detach().numpy()

    #         print("VAL PRED SHAPE", pred.shape)

    #         low_res_dens_path =  os.path.join(
    #                         toy_dir,
    #                         complex_name_val,
    #                         complex_name_val + base_low_res_density_filename
    #                         )
            
    #         low_res_dens_origin = np.array([0, 0, 0])

    #         with mrcfile.open(low_res_dens_path, "r") as low_res_dens_file:
    #             low_res_dens_origin = low_res_dens_file.header.origin

    #         generated_map_path = os.path.join(
    #                         toy_dir,
    #                         complex_name_val,
    #                         complex_name_val + "_val_" + model_name + "_generated_map.mrc"
    #                         )
            
    #         with mrcfile.new(generated_map_path, overwrite=True) as generated_file:
    #             generated_file.set_data(pred)
    #             generated_file.header.origin = low_res_dens_origin
    #             generated_file.voxel_size = (0.5, 0.5, 0.5)


    #     break
