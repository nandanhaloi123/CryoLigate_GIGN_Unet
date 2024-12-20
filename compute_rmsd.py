import numpy as np
import os
import sys
import getopt
from itertools import repeat
from multiprocessing.pool import Pool
from utils import (
    read_molecule,
    find_pdb_ligand_file_in_db,
    find_txt_file_in_db,
    find_file_in_db,
    extract_filename_from_full_path,
    delete_extension_from_filename,
    read_complexes_names_from_file,
    apply_args_and_kwargs,
)

def compute_rmsd(mol1, mol2):

    # Get conformers
    conf1 = mol1.GetConformer()
    conf2 = mol2.GetConformer()

    # Extract coordinates
    coords1 = conf1.GetPositions()
    coords2 = conf2.GetPositions()

    # Compute RMSD
    diff = coords1 - coords2
    rmsd = np.sqrt(np.mean(np.sum(diff**2, axis=1)))

    return rmsd



def main(
        complex_name, 
        db_path, 
        base_ligand_filename="_ligand.pdb",
        base_corr_filename="_corr0.6_passed_complexes.txt",
        ):
    rmsd_list = []
    try:
        corr_path_full = find_txt_file_in_db(complex_name, db_path, base_txt_name=base_corr_filename)
        lignad_path_full = find_pdb_ligand_file_in_db(complex_name, db_path, base_ligand_name=base_ligand_filename)
        ligand_mol = read_molecule(lignad_path_full, remove_Hs=True, sanitize=True)
        with open(corr_path_full, "r") as corr_file:
            for line in corr_file:
                try:
                    rdkit_path_full = line.strip()
                    rdkit_filename = extract_filename_from_full_path(rdkit_path_full)
                    docked_pose_name = delete_extension_from_filename(
                        "_".join(rdkit_filename.split("_")[3:])
                    ) + ".pdb"
                    docked_path_full = find_file_in_db(complex_name, db_path, docked_pose_name)
                    docked_mol = read_molecule(docked_path_full, remove_Hs=True, sanitize=True)
                    rmsd = compute_rmsd(ligand_mol, docked_mol)
                    rmsd_list.append(rmsd)
                except Exception as e:
                    print(f"Failed to compute rmsd for docked pose {docked_pose_name} of complex {complex_name}: {e}")  

    except Exception as e:
        print(f"Failed to compute rmsd for complex {complex_name}: {e}")

    return rmsd_list


    
if __name__ == "__main__":

    # read script's arguments
    opts, args = getopt.getopt(sys.argv[1:], "p:")
    num_process = None  # number of processes for multiprocessing (Pool.starmap)
    for opt, arg in opts:
        if opt == "-p":
            num_process = int(arg)
    assert num_process is not None, "Number of processes is None after arugment's parsing"

    # path to the database with molecule data
    db_path = os.path.sep + os.path.sep.join(
        [
            "mnt",
            "cephfs",
            "projects",
            "2023110101_Ligand_fitting_to_EM_maps",
            "PDBbind",
            "PDBBind_Zenodo_6408497",
        ]
    )

    # load molecule names
    complex_names_csv = (
        db_path + os.path.sep + "PDB_IDs_with_rdkit_length_less_than_16A_succ_gnina.csv"
    )
    complex_names = read_complexes_names_from_file(complex_names_csv)
    n_complexes = len(complex_names)

    base_ligand_filename="_ligand.pdb"
    base_corr_filename="_corr0.6_passed_complexes.txt"

    main_kwargs = {
        "base_ligand_filename": base_ligand_filename,
        "base_corr_filename": base_corr_filename,
    }

    # run computations in several Processes to speed up
    # generate args and kwargs iterables for Pool.starmap()
    main_args_iter = zip(
        complex_names,
        repeat(db_path, n_complexes),
    )  # iterable for positional arguments
    main_kwargs_iter = repeat(main_kwargs, n_complexes)
    starmap_args_iter = zip(
        repeat(main, n_complexes), main_args_iter, main_kwargs_iter
    )  # combined positional and keyword arguments for Pool.starmap()

    rmsd_res_first = []
    rmsd_res_all = [] 
    with Pool(num_process) as pool:
        result = pool.starmap(apply_args_and_kwargs, starmap_args_iter)
        for rmsd_list in result:
            if rmsd_list:
                rmsd_res_first.append(rmsd_list[0])
                for rmsd in rmsd_list:
                    rmsd_res_all.append(rmsd)

    rmsd_res_first = np.array(rmsd_res_first)
    rmsd_res_all = np.array(rmsd_res_all)

    print(f"Average RMSD first: {np.average(rmsd_res_first)}")

    print(f"Average RMSD all: {np.average(rmsd_res_all)}")