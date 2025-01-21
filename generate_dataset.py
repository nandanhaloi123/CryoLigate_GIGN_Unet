# %%
import os
import sys
import getopt
import pandas as pd
import numpy as np
from scipy.spatial import distance_matrix
from itertools import repeat
import torch 
from datetime import datetime, timezone
from multiprocessing.pool import Pool
from torch.utils.data import Dataset, DataLoader
from rdkit import RDLogger
from torch_geometric.data import Batch, Data
import warnings
import mrcfile
RDLogger.DisableLog('rdApp.*')
np.set_printoptions(threshold=np.inf)
warnings.filterwarnings('ignore')

from utils import(
    delete_extension_from_filename,
    log,
    apply_args_and_kwargs,
    find_mrc_density_file_in_db,
    create_folder,
    extract_filename_from_full_path,
    extract_folder_from_full_path,
    find_pdb_ligand_file_in_db,
    find_file_in_db,
)

def normalize_density_minmax(density):
    """
    Normalizes given density array:
    1) Replaces all negative elements with zero
    2) Applies min-max normalization

    Args:
        density - array with density values

    Returns:
        density_normalized - array with normalized density values
    """

    # replace all negative elements with zero
    density[density < 0] = 0

    # apply min-max normalization
    min_value = np.min(density)
    max_value = np.max(density)
    density_normalized = (density - min_value) / (max_value - min_value)

    return density_normalized

def find_top_density(density, percentile):
    """
    Finds top percentile for the given density array.

    Args:
        density - array with density values
        percentile - percentile value between 0 and 1

    Returns:
        density value corresponding the percentile
    """
    use_density = density[density > 0]

    hist, bin_edges = np.histogram(use_density, bins=200)

    log_hist = [np.log(x) if x > 0 else 0 for x in hist]
    sum_cutoff = np.sum(log_hist) * percentile
    cumulative = 0
    for j in range(len(log_hist)):
        cumulative += log_hist[j]
        if cumulative > sum_cutoff:
            return bin_edges[j]
        
    return bin_edges[-1]

def normalize_density_with_percentile(density, percentile=0.99):
    """
    Normalizes given density array:
    1) Replaces all negative elements with zero
    2) Find density value corresponding to the given percentile and set is as a max value
    3) Applies min-max normalization

    Args:
        density - array with density values
        percentile - percentile value between 0 and 1

    Returns:
        density_normalized - array with normalized density values
    """
    # replace all negative elements with zero
    density[density < 0] = 0

    # find percentile density value and set is a max value
    percentile_density_value = find_top_density(density, percentile)
    density[density > percentile_density_value] = percentile_density_value

    # apply min-max normalization
    min_value = np.min(density)
    max_value = np.max(density)
    density_normalized = (density - min_value) / (max_value - min_value)

    return density_normalized


# %%
def generate_network_data(ligand_embedding_path, label_path, low_res_density_path, save_path):
    """
    Generates data for the network training and saves it to the specified file.

    Args:
       ligand_embedding_path - full path to the file with liand embedding (including its name) 
       label_path - full path to the file with target good resolution density (including its name) 
       low_res_density_path - full path to the file with input low resolution density (including its name) 
       save_path - full path to the file where generared data will be saved
    """

    # read ligand embedding data
    ligand_embedding = torch.load(ligand_embedding_path)

    # read and normalize output (target) good resolution density
    label_density = mrcfile.read(label_path)
    label_density_normalized = normalize_density_minmax(label_density)
    y = torch.tensor(label_density_normalized).unsqueeze(0).unsqueeze(0)

    # read and normalize input low resolution density
    low_res_density = mrcfile.read(low_res_density_path)
    low_res_density_normalized = normalize_density_minmax(low_res_density)
    low_res_dens = torch.tensor(low_res_density_normalized).unsqueeze(0).unsqueeze(0)
    
    # save generated data to a file
    data = Data(
        y=y, 
        ligand_embedding=ligand_embedding,
        low_res_dens=low_res_dens
    )

    torch.save(data, save_path)


def get_data(
    complex_name,
    db_path,
    base_ligand_filename="_ligand.pdb",
    base_ligand_embedding_filename="_ligand_embedding.pyg",
    base_label_filename="_ligand_res1.0_boxed_16A.mrc",
    base_low_res_density_filename="_nconfs10_genmodedocking_res3.5_nbox16_threshcorr0.3_delprob0.2_low_resolution_forward_model.mrc",
    clarifying_data_name="",
    create=False,
    is_log=True,
    log_filename="log.txt",
    log_path=os.path.join(os.getcwd(), "dataset_GIGN_main_logs"),
):
    """
    If create==False just checks if the network data exists and returns path to the file.
    If create==True creates data for the network and returns path to the created file.

    Args:
        complex_name - name (id) of the protein-ligand complex
        db_path - path to the database with the complexes' data
        base_ligand_filename - base name for the ligand file (used to construct the full name)
        base_ligand_embedding_filename - base name for the file with the ligand embedding (used to construct the full name)
        base_label_filename - base name for the file with label (target) density (used to construct the full name)
        base_low_res_density_filename - base name for the file with low resolution density (used to construct the full name)
        clarifying_data_name - additional part of the data file name to clarify 
        (needed because we can generate different data files for the same complex e.g. with different good/bad densities)
        create - whether we should create the data file
        is_log - whether we should write logs for the function
        log_filename - name of the log file (needed only if is_log==True)
        log_path - path to the log file excluding its name (needed only if is_log==True)

    Returns:
        full paths to found/created file
    """

    if is_log:
        log(
            f"{"Create" if create else "Look for"} network data for {complex_name}. Db path: {db_path}.",
            status="INFO",
            log_path=log_path,
            log_filename=log_filename,
        )

    if create: # create the data if required
        path_to_created_file = None
        try:
            ligand_path = find_pdb_ligand_file_in_db(
                                    complex_name,
                                    db_path,
                                    base_ligand_name=base_ligand_filename 
                                ) # path to the ligand file
            ligand_embedding_path = find_file_in_db(
                                    complex_name,
                                    db_path,
                                    complex_name + base_ligand_embedding_filename
                                ) # path to the file with ligand embeddings
            label_path = find_mrc_density_file_in_db(
                                    complex_name,
                                    db_path,
                                    base_density_name=base_label_filename 
                                ) # path to the target (good resolution) densities
            low_res_density_path = find_mrc_density_file_in_db(
                                    complex_name,
                                    db_path,
                                    base_density_name=base_low_res_density_filename
                                ) # path to the input (low resolution) densities
            path_to_created_file = create_data(
                            ligand_path,
                            ligand_embedding_path,
                            label_path,
                            low_res_density_path,
                            clarifying_data_name=clarifying_data_name,
                            is_log=is_log,
                            log_filename=log_filename,
                            log_path=log_path,
                            ) 
            
        except Exception as e:
            if is_log:
                log(
                    f"Failed to call create_data() function for the complex {complex_name}: {e}! Db path: {db_path}.",
                    status="ERROR",
                    log_path=log_path,
                    log_filename=log_filename,
                )
        return path_to_created_file

    else: # or try to find data
        path_to_found_file = None
        try:
            ligand_path = find_pdb_ligand_file_in_db(
                        complex_name,
                        db_path,
                        base_ligand_name=base_ligand_filename 
                    ) # path to the ligand file
            path_to_found_file = find_data(
                            ligand_path,
                            clarifying_data_name=clarifying_data_name,
                            is_log=is_log,
                            log_filename=log_filename,
                            log_path=log_path,
                            )
        except Exception as e:
            if is_log:
                log(
                    f"Failed to call find_data() function or the complex {complex_name}: {e}! Db path: {db_path}.",
                    status="ERROR",
                    log_path=log_path,
                    log_filename=log_filename,
                )
        return path_to_found_file 


def create_data(    
    ligand_path,
    ligand_embedding_path,
    label_path,
    low_res_density_path,
    clarifying_data_name="",
    is_log=True,
    log_filename="log.txt",
    log_path=os.path.join(os.getcwd(), "dataset_GIGN_main_logs"),
):
    """
    Creates network data for the complex and returns path to the created file.

    Args:
        ligand_path - full path to the ligand file
        ligand_embedding_path - full path to the file with liand embedding (including its name) 
        label_path - full path to the label (target) good resolution density 
        low_res_density_path - full path to the input low resolution density 
        clarifying_data_name - additional part of the data file name to clarify 
        (needed because we can generate different data files for the same complex e.g. with different good/bad densities)
        is_log - whether we should write logs for the function
        log_filename - name of the log file (needed only if is_log==True)
        log_path - path to the log file excluding its name (needed only if is_log==True)

    Returns:
        path_to_data_file - full path to the network data file or None if failed to create
    """

    ligand_fname = extract_filename_from_full_path(ligand_path) # name of the ligand file
    ligand_folder = extract_folder_from_full_path(ligand_path) # path to the folder with ligand file
    path_to_data_file = os.path.join(
        ligand_folder, 
        f"{delete_extension_from_filename(ligand_fname) + clarifying_data_name}.pyg"
    ) # path to the output file with network data

    # create the network data file
    try:
        generate_network_data(
            ligand_embedding_path,
            label_path,
            low_res_density_path,
            path_to_data_file 
        )
        if is_log:
            log(
                f"Successfully created network data {delete_extension_from_filename(ligand_fname) + clarifying_data_name}.",
                status="INFO",
                log_path=log_path,
                log_filename=log_filename,
            )
        return path_to_data_file
        
    except Exception as e:
        if is_log:
            log(
                f"Failed to create network data {delete_extension_from_filename(ligand_fname) + clarifying_data_name}: {e}",
                status="ERROR",
                log_path=log_path,
                log_filename=log_filename,
            )
        return None


def find_data(   
    ligand_path,  
    clarifying_data_name="",
    is_log=True,
    log_filename="log.txt",
    log_path=os.path.join(os.getcwd(), "dataset_GIGN_main_logs"),
):
    """
    Checks if the data file for the given complex exists and if yes, returns full path to the data file.

    Args:
        ligand_path - full path to the ligand file
        clarifying_data_name - additional part of the data file name to clarify 
        (needed because we can generate different data files for the same complex e.g. with different good/bad densities)
        is_log - whether we should write logs for the function
        log_filename - name of the log file (needed only if is_log==True)
        log_path - path to the log file excluding its name (needed only if is_log==True)

    Returns:
        path_to_data_file - full path to the network data file or None if failed to find
    """
    ligand_fname = extract_filename_from_full_path(ligand_path) # name of the ligand file
    ligand_folder = extract_folder_from_full_path(ligand_path) # path to the folder with ligand file
    path_to_data_file = os.path.join(
        ligand_folder, 
        f"{delete_extension_from_filename(ligand_fname) + clarifying_data_name}.pyg"
    ) # path to the output file with network data

    # look for the data file
    if os.path.isfile(path_to_data_file):
        if is_log:
            log(
                f"Successfully found network data {delete_extension_from_filename(ligand_fname) + clarifying_data_name}.",
                status="INFO",
                log_path=log_path,
                log_filename=log_filename,
            )   
        return path_to_data_file   
     
    else:
        if is_log:
            log(
                f"Failed to find network data {delete_extension_from_filename(ligand_fname) + clarifying_data_name}.",
                status="ERROR",
                log_path=log_path,
                log_filename=log_filename,
        )
        return None


# %%
class PLIDataLoader(DataLoader):
    def __init__(self, data, **kwargs):
        super().__init__(data, collate_fn=data.collate_fn, **kwargs)

class NetworkDataset(Dataset):
    """
    This class is used for generating network data using multiprocessing
    """
    def __init__(
        self, 
        data_dir, 
        data_df, 
        num_process=8, 
        create=False,
        base_ligand_filename="_ligand.pdb",
        base_ligand_embedding_filename="_ligand_embedding.pyg",
        base_label_filename="_ligand_res1.0_boxed_16A.mrc",
        base_low_res_density_filename="_nconfs10_genmodedocking_res3.5_nbox16_threshcorr0.3_delprob0.2_low_resolution_forward_model.mrc",
        clarifying_data_name="",
        is_log=True,
        log_path=os.path.join(os.getcwd(), "dataset_GIGN_main_logs"),
    ):
        self.data_dir = data_dir
        self.data_df = data_df
        self.num_process = num_process
        self.create = create
        self.base_ligand_filename = base_ligand_filename
        self.base_ligand_embedding_filename = base_ligand_embedding_filename
        self.base_label_filename = base_label_filename
        self.base_low_res_density_filename = base_low_res_density_filename
        self.clarifying_data_name = clarifying_data_name
        self.is_log = is_log
        self.log_path = log_path

        self.data_paths = []
        self.complex_ids = []
        self._pre_process()

    def _pre_process(self):
        self.complex_ids = [row['pdbid'] for i, row in self.data_df.iterrows()]
        n_cids = len(self.complex_ids)

        # construct args and kwargs for get_data() function
        data_args_iter = zip(
            self.complex_ids,
            repeat(self.data_dir, n_cids),
        )  # iterable for positional arguments

        log_filename = None
        if self.is_log:
            log_filename = (
            datetime.now(timezone.utc).strftime("%d_%m_%Y_%H.%M.%S.%f")[:-2]
            + "_main_log.txt"
            )
            create_folder(self.log_path)
            
        data_kwargs = {
            "base_ligand_filename": self.base_ligand_filename,
            "base_ligand_embedding_filename": self.base_ligand_embedding_filename,
            "base_label_filename": self.base_label_filename,
            "base_low_res_density_filename": self.base_low_res_density_filename,
            "clarifying_data_name": self.clarifying_data_name,
            "create": self.create,
            "is_log": self.is_log,
            "log_filename": log_filename,
            "log_path": self.log_path
        }
        data_kwargs_iter = repeat(data_kwargs, n_cids) # iterable for keyword arguments
        starmap_args_iter = zip(
            repeat(get_data, n_cids), data_args_iter, data_kwargs_iter
        )  # combined positional and keyword arguments for Pool.starmap()


        # run computations in multiple processes
        with Pool(self.num_process) as pool:
            result = pool.starmap(apply_args_and_kwargs, starmap_args_iter)

        for path_to_data in result:
            if path_to_data is not None:
                self.data_paths.append(path_to_data)


    def __getitem__(self, idx):
        return torch.load(self.data_paths[idx])

    def collate_fn(self, batch):
        return Batch.from_data_list(batch)

    def __len__(self):
        return len(self.data_paths)

if __name__ == '__main__':

    # read script's arguments
    opts, args = getopt.getopt(sys.argv[1:], "p:")
    num_process = None  # number of processes for multiprocessing (Pool.starmap)
    for opt, arg in opts:
        if opt == "-p":
            num_process = int(arg)
    assert num_process is not None, "Number of processes is None after arugment's parsing"

    data_dir = '/proj/berzelius-2022-haloi/users/x_elima/PDBBind_Zenodo_6408497' # path to the database

    # load molecule names
    complex_names_csv = (
        data_dir + os.path.sep + "PDB_IDs_with_rdkit_length_less_than_24A_nbonds_less_than_10.csv"
    )
    data_df = pd.read_csv(complex_names_csv)
    print(f"Generating network data for {len(data_df)} complexes....")

    # generate dataset for the network
    create = True
    base_ligand_filename = "_ligand.pdb"
    base_ligand_embedding_filename = "_ligand_embedding.pyg"
    base_label_filename = "_ligand_res_2.0_gridpsace_0.5_nbox_48_size_24A_label.mrc"
    base_low_res_density_filename = "_nconfs3_genmode_gnina_docking_boxextens1.0_res4.0_delprob0.0_low_resolution_forward_model.mrc"
    clarifying_data_name = "_forward_model_nconfs3_to_good_res2.0_norm_minmax_with_embeddings"
    is_log = True
    log_path = os.path.join(os.getcwd(), "network_dataset_main_logs")
    dataset = NetworkDataset(
        data_dir,
        data_df,
        num_process=num_process,
        create=create,
        base_ligand_filename=base_ligand_filename,
        base_ligand_embedding_filename=base_ligand_embedding_filename,
        base_label_filename=base_label_filename,
        base_low_res_density_filename=base_low_res_density_filename,
        clarifying_data_name=clarifying_data_name,
        is_log=is_log,
        log_path=log_path
    )

    print(len(dataset.data_paths))
    print(len(dataset))
    print(len(np.unique(dataset.data_paths)))
    print(dataset.data_paths)
