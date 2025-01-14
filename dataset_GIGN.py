# %%
import os
import sys
import getopt
import pandas as pd
import numpy as np
import pickle
from scipy.spatial import distance_matrix
import multiprocessing
from itertools import repeat
import networkx as nx
import torch 
from datetime import datetime, timezone
from multiprocessing.pool import Pool
from torch.utils.data import Dataset, DataLoader
from rdkit import Chem
from rdkit import RDLogger
from rdkit import Chem
from torch_geometric.data import Batch, Data
import warnings
import mrcfile
RDLogger.DisableLog('rdApp.*')
np.set_printoptions(threshold=np.inf)
warnings.filterwarnings('ignore')

from utils import(
    find_txt_file_in_db,
    delete_extension_from_filename,
    log,
    apply_args_and_kwargs,
    find_mrc_density_file_in_db,
    create_folder,
)

# %%
def one_of_k_encoding(k, possible_values):
    if k not in possible_values:
        raise ValueError(f"{k} is not a valid value in {possible_values}")
    return [k == e for e in possible_values]


def one_of_k_encoding_unk(x, allowable_set):
    if x not in allowable_set:
        x = allowable_set[-1]
    return list(map(lambda s: x == s, allowable_set))


def atom_features(mol, graph, atom_symbols=['C', 'N', 'O', 'S', 'F', 'P', 'Cl', 'Br', 'I'], explicit_H=True):

    for atom in mol.GetAtoms():
        results = one_of_k_encoding_unk(atom.GetSymbol(), atom_symbols + ['Unknown']) + \
                one_of_k_encoding_unk(atom.GetDegree(),[0, 1, 2, 3, 4, 5, 6]) + \
                one_of_k_encoding_unk(atom.GetImplicitValence(), [0, 1, 2, 3, 4, 5, 6]) + \
                one_of_k_encoding_unk(atom.GetHybridization(), [
                    Chem.rdchem.HybridizationType.SP, Chem.rdchem.HybridizationType.SP2,
                    Chem.rdchem.HybridizationType.SP3, Chem.rdchem.HybridizationType.
                                        SP3D, Chem.rdchem.HybridizationType.SP3D2
                    ]) + [atom.GetIsAromatic()]
        # In case of explicit hydrogen(QM8, QM9), avoid calling `GetTotalNumHs`
        if explicit_H:
            results = results + one_of_k_encoding_unk(atom.GetTotalNumHs(),
                                                    [0, 1, 2, 3, 4])

        atom_feats = np.array(results).astype(np.float32)

        graph.add_node(atom.GetIdx(), feats=torch.from_numpy(atom_feats))

def get_edge_index(mol, graph):
    for bond in mol.GetBonds():
        i = bond.GetBeginAtomIdx()
        j = bond.GetEndAtomIdx()

        graph.add_edge(i, j)

def mol2graph(mol):
    graph = nx.Graph()
    atom_features(mol, graph)
    get_edge_index(mol, graph)

    graph = graph.to_directed()
    x = torch.stack([feats['feats'] for n, feats in graph.nodes(data=True)])
    edge_index = torch.stack([torch.LongTensor((u, v)) for u, v in graph.edges(data=False)]).T

    return x, edge_index

def inter_graph(ligand, pocket, dis_threshold = 5.):
    atom_num_l = ligand.GetNumAtoms()
    atom_num_p = pocket.GetNumAtoms()

    graph_inter = nx.Graph()
    pos_l = ligand.GetConformers()[0].GetPositions()
    pos_p = pocket.GetConformers()[0].GetPositions()
    dis_matrix = distance_matrix(pos_l, pos_p)
    node_idx = np.where(dis_matrix < dis_threshold)
    for i, j in zip(node_idx[0], node_idx[1]):
        graph_inter.add_edge(i, j+atom_num_l) 

    graph_inter = graph_inter.to_directed()
    edge_index_inter = torch.stack([torch.LongTensor((u, v)) for u, v in graph_inter.edges(data=False)]).T

    return edge_index_inter

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
def mols2graphs(complex_path, label_path, low_res_density_path, save_path, dis_threshold=5.):

    with open(complex_path, 'rb') as f:
        ligand, pocket = pickle.load(f)

    atom_num_l = ligand.GetNumAtoms()
    atom_num_p = pocket.GetNumAtoms()

    pos_l = torch.FloatTensor(ligand.GetConformers()[0].GetPositions())
    pos_p = torch.FloatTensor(pocket.GetConformers()[0].GetPositions())
    x_l, edge_index_l = mol2graph(ligand)
    x_p, edge_index_p = mol2graph(pocket)
    x = torch.cat([x_l, x_p], dim=0)
    edge_index_intra = torch.cat([edge_index_l, edge_index_p+atom_num_l], dim=-1)
    edge_index_inter = inter_graph(ligand, pocket, dis_threshold=dis_threshold)

    # read and normalize label (good resolution) density
    # NOTE TODO: change normalizations
    label_density = mrcfile.read(label_path)
    label_density_normalized = normalize_density_minmax(label_density)
    # label_density_normalized = normalize_density_with_percentile(label_density, percentile=0.99)
    y = torch.tensor(label_density_normalized).unsqueeze(0).unsqueeze(0)

    # # NOTE TODO: add normalization back
    # label_density = mrcfile.read(label_path)
    # y = torch.tensor(label_density).unsqueeze(0).unsqueeze(0)

    # read and normalize input (bad resolution) density
    # NOTE TODO: change normalizations
    low_res_density = mrcfile.read(low_res_density_path)
    low_res_density_normalized = normalize_density_minmax(low_res_density)
    # low_res_density_normalized = normalize_density_with_percentile(low_res_density, percentile=0.99)
    low_res_dens = torch.tensor(low_res_density_normalized).unsqueeze(0).unsqueeze(0)

    # # NOTE TODO: add normalization back
    # low_res_density = mrcfile.read(low_res_density_path)
    # low_res_dens = torch.tensor(low_res_density).unsqueeze(0).unsqueeze(0)

    pos = torch.concat([pos_l, pos_p], dim=0)
    split = torch.cat([torch.zeros((atom_num_l, )), torch.ones((atom_num_p,))], dim=0)
    
    data = Data(
        x=x, 
        edge_index_intra=edge_index_intra, 
        edge_index_inter=edge_index_inter, 
        y=y, 
        pos=pos, 
        split=split, 
        low_res_dens=low_res_dens
    )

    torch.save(data, save_path)
    # return data


def get_graphs(
    complex_name,
    db_path,
    base_corr_check_filename="_corr0.6_passed_complexes.txt",
    base_label_filename="_ligand_res1.0_boxed_16A.mrc",
    base_low_res_density_filename="_nconfs10_genmodedocking_res3.5_nbox16_threshcorr0.3_delprob0.2_low_resolution_forward_model.mrc",
    clarifying_graph_name="",
    create=False,
    dis_threshold=5,
    is_log=True,
    log_filename="log.txt",
    log_path=os.path.join(os.getcwd(), "dataset_GIGN_main_logs"),
):
    """
    If create==False just checks if the graphs for complexes from corr_check_file exist and returns those that exist.
    If create==True creates graphs for complexes from corr_check_file and returns paths to those that were successfully created.

    Args:
        complex_name - name (id) of the protein-ligand complex
        db_path - path to the database with the complexes' data
        base_corr_check_filename - base name for the file with complexes that passed cross-correlation check
        (used to construct the full name)
        base_label_filename - base name for the file with label (target) density (used to construct the full name)
        base_low_res_density_filename - base name for the file with low resolution density (used to construct the full name)
        clarifying_graph_name - additional part of the graph name to clarify 
        (needed because we can generate/have different graphs for the same complex e.g. with different good/bad densities)
        create - whether we should create the graphs
        is_log - whether we should write logs for the function
        log_filename - name of the log file (needed only if is_log==True)
        log_path - path to the log file excluding its name (needed only if is_log==True)

    Returns:
        list of full paths to found/created graphs
    """

    if is_log:
        log(
            f"{"Create" if create else "Look for"} graphs for {complex_name}. Db path: {db_path}.",
            status="INFO",
            log_path=log_path,
            log_filename=log_filename,
        )

    if create: # create the graphs if required
        graphs_created = []
        try:
            corr_check_file_path = find_txt_file_in_db(
                                    complex_name,
                                    db_path,
                                    base_txt_name=base_corr_check_filename 
                                ) # path to the complexes passed cross-correlation check
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
            graphs_created = create_graphs(
                            corr_check_file_path,
                            label_path,
                            low_res_density_path,
                            clarifying_graph_name=clarifying_graph_name,
                            dis_threshold=dis_threshold,
                            is_log=is_log,
                            log_filename=log_filename,
                            log_path=log_path,
                            ) 
        except Exception as e:
            if is_log:
                log(
                    f"Failed to call create_graphs() function for the complex {complex_name}: {e}! Db path: {db_path}.",
                    status="ERROR",
                    log_path=log_path,
                    log_filename=log_filename,
                )
        return graphs_created

    else: # or try to find them
        graphs_found = []
        try:
            corr_check_file_path = find_txt_file_in_db(
                                    complex_name,
                                    db_path,
                                    base_txt_name=base_corr_check_filename 
                                ) # path to the complexes passed cross-correlation check
            graphs_found = find_graphs(
                            corr_check_file_path,
                            clarifying_graph_name=clarifying_graph_name,
                            is_log=is_log,
                            log_filename=log_filename,
                            log_path=log_path,
                            )
        except Exception as e:
            if is_log:
                log(
                    f"Failed to call find_graphs() function or the complex {complex_name}: {e}! Db path: {db_path}.",
                    status="ERROR",
                    log_path=log_path,
                    log_filename=log_filename,
                )
        return graphs_found


def create_graphs(    
    corr_check_file_path,
    label_path,
    low_res_density_path,
    clarifying_graph_name="",
    dis_threshold=5,
    is_log=True,
    log_filename="log.txt",
    log_path=os.path.join(os.getcwd(), "dataset_GIGN_main_logs"),
):
    """
    Creates graphs for complexes from corr_check_file and returns paths to those that were successfully created.

    Args:
        corr_check_file_path - full path to the file with path to the complexes that passed cross-correlation check
        label_path - full path to the label (target) density for the complexes
        low_res_density_path - full path to the low resolution density for the complexes
        clarifying_graph_name - additional part of the graph name to clarify 
        (needed because we can generate different graphs for the same complex e.g. with different good/bad densities)
        is_log - whether we should write logs for the function
        log_filename - name of the log file (needed only if is_log==True)
        log_path - path to the log file excluding its name (needed only if is_log==True)
    Returns:
        graphs_created - list with full paths to the created graphs
    """

    paths_to_complexes = [] # list to store paths to the complexes for which graphs should be created
    graphs_to_create = [] # list to store paths to the graphs that should be created
    graphs_created = [] # list to store paths to graphs that were succesfully created 

    # read complexes that passed cross-correlation check
    with open(corr_check_file_path, "r") as complex_file:
        for line in complex_file:
            line = line.strip()
            paths_to_complexes.append(line)
            graphs_to_create.append(delete_extension_from_filename(line) + clarifying_graph_name + ".pyg")

    # create graphs for those complexes
    for i in range(len(paths_to_complexes)):
        try:
            complex_path = paths_to_complexes[i]
            graph_path = graphs_to_create[i]
            mols2graphs(
                complex_path,
                label_path,
                low_res_density_path,
                graph_path,
                dis_threshold=dis_threshold
            )
            graphs_created.append(graph_path)
        except Exception as e:
            if is_log:
                log(
                    f"Failed to create graph for complex {complex_path}: {e}",
                    status="ERROR",
                    log_path=log_path,
                    log_filename=log_filename,
                )
    if is_log:
        log(
            f"Created graphs for {len(graphs_created)} out of {len(graphs_to_create)} complexes from {corr_check_file_path}.",
            status="INFO",
            log_path=log_path,
            log_filename=log_filename,
        )   

    return graphs_created


def find_graphs(    
    corr_check_file_path,
    clarifying_graph_name="",
    is_log=True,
    log_filename="log.txt",
    log_path=os.path.join(os.getcwd(), "dataset_GIGN_main_logs"),
):
    """
    Checks if the graphs for complexes from corr_check_file exist and returns those that exist.

    Args:
        corr_check_file_path - full path to the file with path to the complexes that passed cross-correlation check
        clarifying_graph_name - additional part of the graph name to clarify 
        (needed because we can have different graphs for the same complex e.g. with different good/bad densities)
        is_log - whether we should write logs for the function
        log_filename - name of the log file (needed only if is_log==True)
        log_path - path to the log file excluding its name (needed only if is_log==True)
    Returns:
        graphs_found - list with full paths to the graphs that were found
    """

    graphs_to_find = [] # list to store paths to the graphs that should be found
    graphs_found = [] # list to store paths to graphs that were succesfully found 

    # read complexes that passed cross-correlation check
    with open(corr_check_file_path, "r") as complex_file:
        for line in complex_file:
            line = line.strip()
            graphs_to_find.append(delete_extension_from_filename(line) + clarifying_graph_name + ".pyg")
    #NOTE TODO: remove the thing below, we need all graphs
    graphs_to_find = graphs_to_find[:1]
    # look for graphs for those complexes
    for graph_path in graphs_to_find:
        if os.path.isfile(graph_path):
            graphs_found.append(graph_path)
        else:
            if is_log:
                log(
                    f"Failed to find graph {graph_path}",
                    status="ERROR",
                    log_path=log_path,
                    log_filename=log_filename,
            )
    if is_log:
        log(
            f"Found {len(graphs_found)} out of {len(graphs_to_find)} graphs for complexes from {corr_check_file_path}.",
            status="INFO",
            log_path=log_path,
            log_filename=log_filename,
        )   

    return graphs_found


# %%
class PLIDataLoader(DataLoader):
    def __init__(self, data, **kwargs):
        super().__init__(data, collate_fn=data.collate_fn, **kwargs)

class GraphDataset(Dataset):
    """
    This class is used for generating graph objects using multi process
    """
    def __init__(
        self, 
        data_dir, 
        data_df, 
        dis_threshold=5, 
        graph_type='Graph_GIGN', 
        num_process=8, 
        create=False,
        base_corr_check_filename="_corr0.6_passed_complexes.txt",
        base_label_filename="_ligand_res1.0_boxed_16A.mrc",
        base_low_res_density_filename="_nconfs10_genmodedocking_res3.5_nbox16_threshcorr0.3_delprob0.2_low_resolution_forward_model.mrc",
        clarifying_graph_name="",
        is_log=True,
        log_path=os.path.join(os.getcwd(), "dataset_GIGN_main_logs"),
    ):
        self.data_dir = data_dir
        self.data_df = data_df
        self.dis_threshold = dis_threshold
        self.graph_type = graph_type
        self.create = create
        self.graph_paths = []
        self.complex_ids = []
        self.num_process = num_process
        self.base_corr_check_filename = base_corr_check_filename
        self.base_label_filename = base_label_filename
        self.base_low_res_density_filename = base_low_res_density_filename
        self.clarifying_graph_name = clarifying_graph_name
        self.is_log = is_log
        self.log_path = log_path
        self._pre_process()

    def _pre_process(self):
        # dis_thresholds = repeat(self.dis_threshold, len(data_df))

        self.complex_ids = [row['pdbid'] for i, row in self.data_df.iterrows()]
        n_cids = len(self.complex_ids)

        # construct args and kwargs for get_graphs() function
        graph_args_iter = zip(
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
            
        graph_kwargs = {
            "base_corr_check_filename": self.base_corr_check_filename,
            "base_label_filename": self.base_label_filename,
            "base_low_res_density_filename": self.base_low_res_density_filename,
            "clarifying_graph_name": self.clarifying_graph_name,
            "create": self.create,
            "dis_threshold": self.dis_threshold,
            "is_log": self.is_log,
            "log_filename": log_filename,
            "log_path": self.log_path
        }
        graph_kwargs_iter = repeat(graph_kwargs, n_cids) # iterable for keyword arguments
        starmap_args_iter = zip(
            repeat(get_graphs, n_cids), graph_args_iter, graph_kwargs_iter
        )  # combined positional and keyword arguments for Pool.starmap()

        # run computations in multiple processes
        with Pool(self.num_process) as pool:
            result = pool.starmap(apply_args_and_kwargs, starmap_args_iter)

        for res in result:
            for path in res:
                self.graph_paths.append(path)


    def __getitem__(self, idx):
        return torch.load(self.graph_paths[idx])

    def collate_fn(self, batch):
        return Batch.from_data_list(batch)

    def __len__(self):
        return len(self.graph_paths)

if __name__ == '__main__':
    # read script's arguments
    opts, args = getopt.getopt(sys.argv[1:], "p:")
    num_process = None  # number of processes for multiprocessing (Pool.starmap)
    for opt, arg in opts:
        if opt == "-p":
            num_process = int(arg)
    assert num_process is not None, "Number of processes is None after arugment's parsing"

    # path to the database
    # data_dir = os.path.sep + os.path.sep.join(
    #     [
    #         "mnt",
    #         "cephfs",
    #         "projects",
    #         "2023110101_Ligand_fitting_to_EM_maps",
    #         "PDBbind",
    #         "PDBBind_Zenodo_6408497",
    #     ]
    # )

    data_dir = '/proj/berzelius-2022-haloi/users/x_elima/PDBBind_Zenodo_6408497'

    # load molecule names
    complex_names_csv = (
        data_dir + os.path.sep + "PDB_IDs_with_rdkit_length_less_than_24A.csv"
    )
    data_df = pd.read_csv(complex_names_csv)
    # data_df = data_df.iloc[:5]
    print(f"Computing graphs for {len(data_df)} complexes....")
    # generate dataset of graphs
    dis_threshold = 5
    graph_type = 'Graph_GIGN' 
    create = True
    base_corr_check_filename = "_complexes_temporary_bad_maps_only.txt"
    base_label_filename = "_ligand_res_2.0_gridpsace_0.5_nbox_48_size_24A_label.mrc"
    base_low_res_density_filename = "_nconfs3_genmode_gnina_docking_boxextens1.0_res4.0_delprob0.0_low_resolution_forward_model.mrc"
    clarifying_graph_name = "_forward_model_bad_nconfs3_to_good_res2.0_norm_minmax"
    is_log = True
    log_path = os.path.join(os.getcwd(), "dataset_GIGN_main_logs")
    dataset = GraphDataset(
        data_dir,
        data_df,
        dis_threshold=dis_threshold,
        graph_type=graph_type,
        num_process=num_process,
        create=create,
        base_corr_check_filename=base_corr_check_filename,
        base_label_filename=base_label_filename,
        base_low_res_density_filename=base_low_res_density_filename,
        clarifying_graph_name=clarifying_graph_name,
        is_log=is_log,
        log_path=log_path
    )

    print(len(dataset.graph_paths))
    print(len(dataset))
    print(len(np.unique(dataset.graph_paths)))
    print(dataset.graph_paths)
