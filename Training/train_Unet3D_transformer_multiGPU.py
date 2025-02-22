import os
import sys
import torch
import torch.nn as nn
import torch.optim as optim
import torch.distributed as dist
import torch.multiprocessing as mp
from torch.nn.parallel import DistributedDataParallel as DDP
import joblib
import pandas as pd
import numpy as np
from datetime import datetime
from sklearn.model_selection import KFold

# Append repo path for convenient imports
sys.path.append(os.path.dirname(os.path.dirname(os.path.realpath(__file__))))

from utils import AverageMeter
from model.Unet3D_transformer import CryoLigate
from data_generation.generate_dataset import NetworkDataset, PLIDataLoader
from config.config_dict import Config
from log.train_logger import TrainLogger

class CustomLoss(nn.Module):
    """Combined loss: MSE + SSIM"""
    def __init__(self):
        super(CustomLoss, self).__init__()

    def forward(self, inputs, targets):
        mse_loss = ((inputs - targets) ** 2).mean()
        epsilon = 1e-6
        n_points = inputs[0].numel()
        inputs_mean = inputs.mean(dim=(2, 3, 4), keepdim=True)
        targets_mean = targets.mean(dim=(2, 3, 4), keepdim=True)
        cov = ((inputs - inputs_mean) * (targets - targets_mean)).sum(dim=(1, 2, 3, 4)) / (n_points - 1)
        inputs_var = inputs.var(dim=(1, 2, 3, 4))
        targets_var = targets.var(dim=(1, 2, 3, 4))
        ssim_loss = (1.0 - (2 * cov + epsilon) / (inputs_var + targets_var + epsilon)).mean()
        return mse_loss + ssim_loss

def val(model, dataloader, device, criterion):
    """Validation step"""
    model.eval()
    val_mse_loss = AverageMeter()
    val_ssim_loss = AverageMeter()
    with torch.no_grad():
        for data in dataloader:
            data = data.to(device)
            pred = model(data)
            label = data.y
            val_mse, val_ssim = criterion.separate_losses(pred, label)
            val_mse_loss.update(val_mse, label.size(0))
            val_ssim_loss.update(val_ssim, label.size(0))
    return val_mse_loss.get_average(), val_ssim_loss.get_average()

def setup(rank, world_size):
    """Initialize DDP environment"""
    os.environ["MASTER_ADDR"] = "localhost"
    os.environ["MASTER_PORT"] = "12355"
    dist.init_process_group("nccl", rank=rank, world_size=world_size)
    torch.cuda.set_device(rank)

def cleanup():
    """Cleanup DDP"""
    dist.destroy_process_group()

def train(rank, world_size):
    """Main training loop"""
    setup(rank, world_size)
    
    # Config
    cfg = 'TrainConfig_CryoLigate'
    config = Config(cfg)
    args = config.get_config()
    batch_size = 256 // world_size  # Adjust batch size per GPU
    epochs = 500
    lr = 5e-4
    wd = 1e-4

    # Load dataset
    data_root = '/proj/berzelius-2022-haloi/users/x_nanha'
    toy_dir = os.path.join(data_root, 'PDBBind_Zenodo_6408497')
    toy_df = pd.read_csv(os.path.join(toy_dir, "PDB_IDs_with_rdkit_length_less_than_24A_nbonds_less_than_10.csv"))
    
    kfold = KFold(n_splits=10, shuffle=True, random_state=42)
    train_idx, val_idx = list(kfold.split(toy_df))[0]
    train_df, valid_df = toy_df.iloc[train_idx], toy_df.iloc[val_idx]

    # DDP DataLoader
    train_set = NetworkDataset(toy_dir, train_df)
    valid_set = NetworkDataset(toy_dir, valid_df)
    train_sampler = torch.utils.data.distributed.DistributedSampler(train_set, num_replicas=world_size, rank=rank)
    train_loader = PLIDataLoader(train_set, batch_size=batch_size, sampler=train_sampler)
    valid_loader = PLIDataLoader(valid_set, batch_size=batch_size)

    # Model
    device = torch.device(f"cuda:{rank}")
    model = CryoLigate().to(device)
    model = DDP(model, device_ids=[rank])

    optimizer = optim.Adam(model.parameters(), lr=lr, weight_decay=wd)
    criterion = CustomLoss()
    best_val_loss = float('inf')

    # Training Loop
    for epoch in range(epochs):
        train_sampler.set_epoch(epoch)
        model.train()
        train_loss_mse, train_loss_ssim = AverageMeter(), AverageMeter()

        for data in train_loader:
            data = data.to(device)
            pred = model(data)
            label = data.y

            loss = criterion(pred, label)
            optimizer.zero_grad()
            loss.backward()
            optimizer.step()

            # Compute train loss
            train_mse, train_ssim = criterion.separate_losses(pred, label)
            train_loss_mse.update(train_mse, label.size(0))
            train_loss_ssim.update(train_ssim, label.size(0))

        # Validation
        val_mse, val_ssim = val(model, valid_loader, device, criterion)
        total_val_loss = val_mse + val_ssim

        # Save best model
        if rank == 0 and total_val_loss < best_val_loss:
            best_val_loss = total_val_loss
            model_dir = os.path.join(os.getcwd(), "saved_models")
            os.makedirs(model_dir, exist_ok=True)
            model_pkl_file = os.path.join(model_dir, f"model_{epoch}.pkl")
            joblib.dump(model.module, model_pkl_file)

        print(f"[Rank {rank}] Epoch {epoch} | Train Loss: {train_loss_mse.get_average():.6f}, {train_loss_ssim.get_average():.6f} | Val Loss: {val_mse:.6f}, {val_ssim:.6f}")

    cleanup()

if __name__ == "__main__":
    world_size = torch.cuda.device_count()
    mp.spawn(train, args=(world_size,), nprocs=world_size, join=True)
