
# %%
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


# %%
# def val(model, dataloader, device):
#     model.eval()

#     pred_list = []
#     label_list = []
#     for data in dataloader:
#         data = data.to(device)
#         with torch.no_grad():
#             pred = model(data)
#             label = data.y

#             pred_list.append(pred.detach().cpu().numpy())
#             label_list.append(label.detach().cpu().numpy())
            
#     pred = np.concatenate(pred_list, axis=0)
#     label = np.concatenate(label_list, axis=0)

#     coff = np.corrcoef(pred, label)[0, 1]
#     # epoch_rmse = criterion(pred, label)
#     epoch_rmse = np.sqrt(mean_squared_error(label, pred))

#     model.train()

#     return epoch_rmse, coff

def val(model, dataloader, device):
    model.eval()

    criterion = nn.MSELoss()
    val_loss = AverageMeter()
    for data in dataloader:
        data = data.to(device)
        with torch.no_grad():
            pred = model(data)
            label = data.y
            loss = criterion(pred, label)
            val_loss.update(loss.item(), label.size(0))

    epoch_loss = val_loss.get_average()
    val_loss.reset()
    model.train()

    return epoch_loss 

# %%
if __name__ == '__main__':
    cfg = 'TrainConfig_GIGN'
    config = Config(cfg)
    args = config.get_config()
    graph_type = args.get("graph_type")
    save_model = args.get("save_model")
    batch_size = 16
    data_root = '/mnt/cephfs/projects/2023110101_Ligand_fitting_to_EM_maps/PDBbind'
    epochs = 300
    repeats = args.get('repeat')
    repeats = 1
    early_stop_epoch = args.get("early_stop_epoch")

    toy_dir = os.path.join(data_root, 'PDBBind_Zenodo_6408497')


    ######################
    toy_df = pd.read_csv(os.path.join(toy_dir, "PDB_IDs_with_rdkit_length_less_than_16A_succ_gnina.csv")).sample(frac=1., random_state=123)
    base_model_name = "only_final_unet_without_final_ReLU_unnorm_maps_new_approach_less_noisy_bad_forward_to_2.0_batchsize_16_hidsize_256_levels_256_lr_5e-4_wd_1e-6"

    split_idx = int(0.9 * len(toy_df))
    train_df = toy_df.iloc[:split_idx]
    valid_df = toy_df.iloc[split_idx:]

    dis_threshold = 5
    num_process = 20
    create_dataset = False
    # base_corr_check_filename = "_Nandan_graphs_with_bad_maps.txt"
    # base_corr_check_filename = "_corr0.6_passed_complexes.txt"
    base_corr_check_filename = "_complexes_map_less_noisy_bad_to_map2.0.txt"
    base_label_filename = "_ligand_res_2.0_gridpsace_0.5_nbox_32_size_16A.mrc"
    base_low_res_density_filename = "_nconfs3_genmode_gnina_docking_boxextens1.0_res4.0_delprob0.0_low_resolution_forward_model.mrc"
    clarifying_graph_name = "_unnormalized_map_new_approach_less_noisy_bad_to_map2.0"
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

    train_loader = PLIDataLoader(train_set, batch_size=batch_size, shuffle=True, drop_last=True)
    valid_loader = PLIDataLoader(valid_set, batch_size=batch_size, shuffle=False, num_workers=4)
    
    model_name = base_model_name
    args["model_name"] = model_name
    logger = TrainLogger(args, cfg, create=True)
    logger.info(__file__)
    logger.info(f"train data: {len(train_set)}")
    logger.info(f"train unqiue data: {len(np.unique(train_set.graph_paths))}")
    logger.info(f"valid data: {len(valid_set)}")
    logger.info(f"valid unqiue data: {len(np.unique(valid_set.graph_paths))}")

    device = torch.device('cuda:0')
    # model = GIGN(35, 512).to(device)
    model = GIGN(35, 256).to(device)
    
    optimizer = optim.Adam(model.parameters(), lr=5e-4, weight_decay=1e-6)
    criterion = nn.MSELoss()

    running_loss = AverageMeter()
    running_acc = AverageMeter()
    running_best_mse = BestMeter("min")
    best_model_list = []

    # initial_train_loss = val(model, train_loader, device)
    # initial_val_loss = val(model, valid_loader, device)

    # msg = "Initial_train_loss-%.7f  Initial_val_loss-%.7f" \
    # % (initial_train_loss, initial_val_loss)
    # logger.info(msg)
    
    best_val_loss = 10000000

    model.train()
    for epoch in range(epochs):
        for data in train_loader:
            data = data.to(device)
            print("Data size", data.size())
            pred = model(data)
            label = data.y
            print("PRED", pred.size())
            print("LABEL", label.size())

            loss = criterion(pred, label)
            optimizer.zero_grad()
            loss.backward()
            optimizer.step()

            running_loss.update(loss.item(), label.size(0)) 

        epoch_train_loss = running_loss.get_average()
        running_loss.reset()
        
        epoch_val_loss = val(model, valid_loader, device)

        msg = "epoch-%d, train_loss-%.7f, val_loss-%.7f" \
        % (epoch, epoch_train_loss, epoch_val_loss)
        logger.info(msg)

        # if (epoch >= 10) and (epoch_val_loss < best_val_loss):
        #     model_dir = logger.get_model_dir()
        #     model_pkl_file = os.path.join(model_dir, "model.pkl")
        #     joblib.dump(model, model_pkl_file)
        #     best_val_loss = epoch_val_loss

        #     msg = "Saved new best model at epoch %d with validation loss %.7f" \
        #     % (epoch, best_val_loss)
        #     logger.info(msg)

        if epoch % 10 == 0:
            model_dir = logger.get_model_dir()
            model_pkl_file = os.path.join(model_dir, f"model_{epoch}.pkl")
            joblib.dump(model, model_pkl_file)

            msg = "Saved model at epoch %d" % (epoch)
            logger.info(msg)


    ##################################
    # toy_df = pd.read_csv(os.path.join(toy_dir, "PDB_IDs_with_rdkit_length_less_than_16A_succ_gnina.csv"))
    # k_folds = 5
    # kfold = KFold(n_splits=k_folds, shuffle=True, random_state=42)
    # for fold, (train_idx, val_idx) in enumerate(kfold.split(toy_df)):
    #     train_df = toy_df.iloc[train_idx]
    #     valid_df = toy_df.iloc[val_idx]
    #     dis_threshold = 5
    #     num_process = 20
    #     create_dataset = False
    #     base_corr_check_filename = "_Nandan_graphs_with_bad_maps.txt"
    #     # base_corr_check_filename = "_Nandan_graphs.txt"
    #     # base_corr_check_filename = "_corr0.6_passed_complexes.txt"
    #     base_label_filename = "_ligand_res_2.0_gridpsace_0.5_nbox_32_size_16A.mrc"
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
    #                 clarifying_graph_name=clarifying_graph_name,
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
    #                 clarifying_graph_name=clarifying_graph_name,
    #                 is_log=is_dataset_log,
    #                 log_path=dataset_log_path,
    #             )
    
    #     train_loader = PLIDataLoader(train_set, batch_size=batch_size, shuffle=False, drop_last=True)
    #     valid_loader = PLIDataLoader(valid_set, batch_size=batch_size, shuffle=False, num_workers=4)
        
    #     model_name = f"k{fold}_" + base_model_name
    #     args["model_name"] = model_name
    #     logger = TrainLogger(args, cfg, create=True)
    #     logger.info(__file__)
    #     logger.info(f"train data: {len(train_set)}")
    #     logger.info(f"valid data: {len(valid_set)}")

    #     device = torch.device('cuda:0')
    #     model = GIGN(35, 512).to(device)
    #     # model = GIGN(35, 256).to(device)
        
    #     optimizer = optim.Adam(model.parameters(), lr=5e-4, weight_decay=1e-6)
    #     criterion = nn.MSELoss()

    #     running_loss = AverageMeter()
    #     running_acc = AverageMeter()
    #     running_best_mse = BestMeter("min")
    #     best_model_list = []

    #     initial_train_loss = val(model, train_loader, device)
    #     initial_val_loss = val(model, valid_loader, device)

    #     msg = "Initial_train_loss-%.7f  Initial_val_loss-%.7f" \
    #     % (initial_train_loss, initial_val_loss)
    #     logger.info(msg)
    #
    #
    #     best_val_loss = 10000000
    #     model.train()
    #     for epoch in range(epochs):
    #         for data in train_loader:
    #             data = data.to(device)
    #             print("Data size", data.size())
    #             pred = model(data)
    #             label = data.y
    #             print("PRED", pred.size())
    #             print("LABEL", label.size())

    #             loss = criterion(pred, label)
    #             optimizer.zero_grad()
    #             loss.backward()
    #             optimizer.step()

    #             running_loss.update(loss.item(), label.size(0)) 

    #         epoch_train_loss = running_loss.get_average()
    #         running_loss.reset()
            
    #         epoch_val_loss = val(model, valid_loader, device)

    #         msg = "epoch-%d, train_loss-%.7f, val_loss-%.7f" \
    #         % (epoch, epoch_train_loss, epoch_val_loss)
    #         logger.info(msg)

            # if (epoch >= 10) and (epoch_val_loss < best_val_loss):
            #     model_dir = logger.get_model_dir()
            #     model_pkl_file = os.path.join(model_dir, "model.pkl")
            #     joblib.dump(model, model_pkl_file)
            #     best_val_loss = epoch_val_loss
                
            #     msg = "Saved new best model at epoch %d with validation loss %.7f" \
            #     % (epoch, best_val_loss)
            #     logger.info(msg)



            # f.write(f"Epoch {epoch}: training loss {epoch_loss}  val loss: {epoch_val_loss}\n")
            # epoch_val_loss = val(model, valid_loader, device)
            # f.write(f"Epoch {epoch}: training loss {epoch_loss}  val loss: {epoch_val_loss}\n")


            # msg = "epoch-%d, train_loss-%.4f, valid_loss-%.4f" \
            # % (epoch, epoch_loss, epoch_val_loss)
            # logger.info(msg)

            # #start validating
            # valid_rmse, valid_pr = val(model, valid_loader, device)
            # msg = "epoch-%d, train_loss-%.4f, train_rmse-%.4f, valid_rmse-%.4f, valid_pr-%.4f" \
            #         % (epoch, epoch_loss, epoch_rmse, valid_rmse, valid_pr)
            # logger.info(msg)

            # if valid_rmse < running_best_mse.get_best():
            #     running_best_mse.update(valid_rmse)
            #     if save_model:
            #         msg = "epoch-%d, train_loss-%.4f, train_rmse-%.4f, valid_rmse-%.4f, valid_pr-%.4f" \
            #         % (epoch, epoch_loss, epoch_rmse, valid_rmse, valid_pr)
            #         model_path = os.path.join(logger.get_model_dir(), msg + '.pt')
            #         best_model_list.append(model_path)
            #         save_model_dict(model, logger.get_model_dir(), msg)
            # else:
            #     count = running_best_mse.counter()
            #     if count > early_stop_epoch:
            #         best_mse = running_best_mse.get_best()
            #         msg = "best_rmse: %.4f" % best_mse
            #         logger.info(f"early stop in epoch {epoch}")
            #         logger.info(msg)
            #         break_flag = True
            #         break


# %%    
