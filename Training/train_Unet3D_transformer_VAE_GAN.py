
# %%
import os
import sys
import joblib
#os.environ['CUDA_VISIBLE_DEVICES'] = "0"
import torch
torch.cuda.empty_cache()
import torch.nn as nn
import torch.optim as optim
import pandas as pd
import numpy as np
from datetime import datetime, timezone
from sklearn.metrics import mean_squared_error
from sklearn.model_selection import KFold


# append repo path to sys for convenient imports
sys.path.append(os.path.dirname(os.path.dirname(os.path.realpath(__file__))))

from utils import AverageMeter
from model.Unet3D_transformer_VAE_GAN import CryoLigateVAE
from model.Unet3D_transformer_VAE_GAN import Discriminator3D
from data_generation.generate_dataset import NetworkDataset, PLIDataLoader
from config.config_dict import Config
from log.train_logger import TrainLogger
from utils import *


class CustomLoss(nn.Module):
    def __init__(self):
        super(CustomLoss, self).__init__()

    def forward(self, inputs, targets, mu, logvar):
        """
        Computes the combined loss: reconstruction loss + KL divergence
        """
        # Reconstruction loss (MSE)
        mse_loss = ((inputs - targets) ** 2).mean()

        # KL Divergence (for VAE)
        # Compute the KL divergence using mu and logvar
        kl_loss = -0.5 * torch.mean(1 + logvar - mu.pow(2) - logvar.exp())

        # Total loss (MSE + KL divergence)
        return mse_loss + kl_loss
    
    def separate_losses(self, inputs, targets, mu, logvar):
        """
        Computes individual losses: MSE and KL divergence
        """
        # Compute MSE loss
        mse_loss = ((inputs - targets) ** 2).mean()

        # Compute KL divergence
        kl_loss = -0.5 * torch.mean(1 + logvar - mu.pow(2) - logvar.exp())
        
        return mse_loss.item(), kl_loss.item()


def val(model, dataloader, device, criterion):
    model.eval()

    val_mse_loss = AverageMeter()
    val_kl_loss = AverageMeter()
    for data in dataloader:
        data = data.to(device)
        with torch.no_grad():
            pred, mu, logvar = model(data)  # Get predictions, mu, and logvar
            label = data.y
            val_mse, val_kl = criterion.separate_losses(pred, label, mu, logvar)
            val_mse_loss.update(val_mse, label.size(0))
            val_kl_loss.update(val_kl, label.size(0))

    # Get average loss values
    epoch_mse_loss = val_mse_loss.get_average()
    epoch_kl_loss = val_kl_loss.get_average()
    val_mse_loss.reset()
    val_kl_loss.reset()

    model.train()

    return epoch_mse_loss, epoch_kl_loss


if __name__ == '__main__':

    # config for the proper model saving
    cfg = 'TrainConfig_CryoLigate'
    config = Config(cfg)
    args = config.get_config()
    batch_size = 16
    epochs = 500
    lr = 5e-4
    wd = 1e-4
    lamd = 0.1
    # data paths
    data_root = '/proj/berzelius-2022-haloi/users/x_nanha'
    toy_dir = os.path.join(data_root, 'PDBBind_Zenodo_6408497')
    toy_df = pd.read_csv(os.path.join(toy_dir, "PDB_IDs_with_rdkit_length_less_than_24A_nbonds_less_than_10.csv"))

    # cross-validation splitting
    n_splits = 10
    kfold = KFold(n_splits=n_splits, shuffle=True, random_state=42)
    kfold_split = list(kfold.split(toy_df))
    k = 0
    train_idx, val_idx = kfold_split[k]
    train_df = toy_df.iloc[train_idx]
    valid_df = toy_df.iloc[val_idx]

    # clear name for the model (to distinguish in the future)
    model_name = f"k_{k}_Unet3D_transformer_VAE_GAN_{lamd}_with_Ligand_embeddings_Hybrid_loss_Norm_minmax_maps_Forward_model_bad_nconfs3_to_Good_res2.0_Batchsize_{batch_size}_lr_{lr:.1e}_wd_{wd:.1e}"
    args["model_name"] = model_name

    # find and read training and validation data
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
    train_loader = PLIDataLoader(train_set, batch_size=batch_size, shuffle=False, drop_last=True)
    valid_loader = PLIDataLoader(valid_set, batch_size=batch_size, shuffle=False, num_workers=0)
    
    # logging some initial information
    logger = TrainLogger(args, cfg, create=True)
    # logger.info(__file__)
    logger.info(f"train data: {len(train_set)}")
    logger.info(f"train unqiue data: {len(np.unique(train_set.data_paths))}")
    logger.info(f"valid data: {len(valid_set)}")
    logger.info(f"valid unqiue data: {len(np.unique(valid_set.data_paths))}")

    # specify the model and optimizer
    # device = torch.device('cuda:0')
    device = torch.device("cuda:0" if torch.cuda.is_available() else "cpu")
    model = CryoLigateVAE().to(device)
    discriminator = Discriminator3D().to(device)

    # if torch.cuda.device_count() > 1:
    #     model = nn.DataParallel(model)
    #     discriminator = nn.DataParallel(discriminator)

    optimizer_G = optim.Adam(model.parameters(), lr=lr, weight_decay=wd)
    optimizer_D = optim.Adam(discriminator.parameters(), lr=lr, weight_decay=wd)

    criterion = CustomLoss()
    adversarial_loss = nn.BCELoss()  # Binary Cross-Entropy for GAN

    val_criterion = CustomLoss()
    train_loss_mse = AverageMeter()
    train_loss_kl = AverageMeter()

    # initial_train_loss = val(model, train_loader, device, val_criterion)
    # initial_val_loss = val(model, valid_loader, device, val_criterion)

    # msg = "Initial_train_loss-%.7f  Initial_val_loss-%.7f" \
    # % (initial_train_loss, initial_val_loss)
    # logger.info(msg)
    
    best_val_loss = 10000000 # initial value of the best validation loss (should be high for the proper model saving)

    model.train()
    discriminator.train()
    # ===== Train Loop =====
for epoch in range(epochs):
    for data in train_loader:
        data = data.to(device)  # Move input data to the same device for both models
        pred, mu, logvar = model(data)  # Generate predictions
        label = data.y.to(device)  # Move label to the same device for discriminator

        # ===== Train Discriminator =====
        optimizer_D.zero_grad()

        real_labels = torch.ones(label.size(0), 1, device=device)
        fake_labels = torch.zeros(label.size(0), 1, device=device)

        real_preds = discriminator(label)  # Discriminator on real data
        fake_preds = discriminator(pred.detach())  # Discriminator on fake data (no gradient for generator)

        loss_real = adversarial_loss(real_preds, real_labels)
        loss_fake = adversarial_loss(fake_preds, fake_labels)
        loss_D = (loss_real + loss_fake) / 2
        loss_D.backward()
        optimizer_D.step()

        # ===== Train Generator (VAE) =====
        optimizer_G.zero_grad()

        vae_loss = criterion(pred, label, mu, logvar)
        fake_preds = discriminator(pred)  # Get discriminator's opinion on generated data

        loss_G = vae_loss + lamd * adversarial_loss(fake_preds, real_labels)  # Adversarial loss encourages realism
        loss_G.backward()
        optimizer_G.step()

    # Validation
    epoch_val_mse, epoch_val_kl = val(model, valid_loader, device, criterion)  # Ensure val uses the same device
    total_val_loss = epoch_val_mse + epoch_val_kl

    # Logging
    msg = f"Epoch-{epoch}, Loss_G: {loss_G.item():.7f}, Loss_D: {loss_D.item():.7f}, Val_MSE: {epoch_val_mse:.7f}, Val_KL: {epoch_val_kl:.7f}"
    logger.info(msg)

    # Save Best Model
    if total_val_loss < best_val_loss:
        model_dir = logger.get_model_dir()
        model_pkl_file = os.path.join(model_dir, f"model_{epoch}.pkl")
        joblib.dump(model, model_pkl_file)
        best_val_loss = total_val_loss
        logger.info(f"Saved new best model at epoch {epoch} with validation loss {best_val_loss:.7f}")
