# %%
import os
import sys
import pickle
import getopt
from itertools import repeat
from datetime import datetime, timezone
from multiprocessing.pool import Pool
from rdkit import Chem
from rdkit.Chem import AllChem
import pandas as pd
from tqdm import tqdm
import pymol
from rdkit import RDLogger
RDLogger.DisableLog('rdApp.*')
from openbabel import openbabel
import numpy as np
import mrcfile

# append repo path to sys for convenient imports
sys.path.append(os.path.dirname(os.path.dirname(os.path.realpath(__file__))))

from utils import (
    delete_extension_from_filename, 
    extract_format_from_filename, 
    extract_filename_from_full_path,
    find_pdb_ligand_file_in_db,
    find_pdb_protein_file_in_db,
    log,
    read_complexes_names_from_file,
    create_folder,
    apply_args_and_kwargs,
)

# %%
# def generate_pocket(data_dir, data_df, distance=5):
#     # complex_id = os.listdir(data_dir)
#     # for cid in complex_id:
#     for i, row in data_df.iterrows():
#         cid = row['pdbid']
#         complex_dir = os.path.join(data_dir, cid)
#         lig_native_path = os.path.join(complex_dir, f"{cid}_ligand.mol2")
#         #protein_path= os.path.join(complex_dir, f"{cid}_protein.pdb")
#         protein_path= os.path.join(complex_dir, f"{cid}_protein_processed.pdb")

#         if os.path.exists(os.path.join(complex_dir, f'Pocket_{distance}A.pdb')):
#             continue

#         pymol.cmd.load(protein_path)
#         pymol.cmd.remove('resn HOH')
#         pymol.cmd.load(lig_native_path)
#         pymol.cmd.remove('hydrogens')
#         pymol.cmd.select('Pocket', f'byres {cid}_ligand around {distance}')
#         pymol.cmd.save(os.path.join(complex_dir, f'Pocket_{distance}A.pdb'), 'Pocket')
#         pymol.cmd.delete('all')

# def generate_complex(data_dir, data_df, distance=5, input_ligand_format='mol2'):
#     pbar = tqdm(total=len(data_df))
#     for i, row in data_df.iterrows():
#         # cid, pKa = row['pdbid'], float(row['-logKd/Ki'])
#         cid = row['pdbid']
#         complex_dir = os.path.join(data_dir, cid)
#         pocket_path = os.path.join(data_dir, cid, f'Pocket_{distance}A.pdb')
#         if input_ligand_format != 'pdb':
#             ligand_input_path = os.path.join(data_dir, cid, f'{cid}_ligand.{input_ligand_format}')
#             ligand_path = ligand_input_path.replace(f".{input_ligand_format}", ".pdb")
#             # os.system(f'obabel {ligand_input_path} -O {ligand_path} -d')
#             obConversion = openbabel.OBConversion()
#             obConversion.SetInAndOutFormats(input_ligand_format, "pdb")
#             mol = openbabel.OBMol()
#             obConversion.ReadFile(mol, ligand_input_path)
#             obConversion.WriteFile(mol, ligand_path)
#         else:
#             ligand_path = os.path.join(data_dir, cid, f'{cid}_ligand.pdb')

#         save_path = os.path.join(complex_dir, f"{cid}_{distance}A.rdkit")
#         ligand = Chem.MolFromPDBFile(ligand_path, removeHs=True)
#         if ligand == None:
#             print(f"Unable to process ligand of {cid}")
#             continue

#         pocket = Chem.MolFromPDBFile(pocket_path, removeHs=True)
#         if pocket == None:
#             print(f"Unable to process protein of {cid}")
#             continue

#         complex = (ligand, pocket)
#         with open(save_path, 'wb') as f:
#             pickle.dump(complex, f)

#         pbar.update(1)



def generate_pocket(ligand_path_full, protein_path_full, output_pocket_path_full, distance=5):
    """
    Generates a pocket for the given ligand and protein using Pymol and saves it to a file.

    Args: 
        ligand_path_full - full path to the input ligand file (including its name)
        protein_path_full - full path to the input protein file (including its name)
        output_pocket_path_full - full path to the output pocket file (including its name)
        distance - threshold distance to generate the pocket
    """
    # get name of the ligand file
    ligand_filename = extract_filename_from_full_path(ligand_path_full)

    # generate pocket using pymol
    pymol.cmd.load(protein_path_full)
    pymol.cmd.remove('resn HOH')
    pymol.cmd.load(ligand_path_full)
    pymol.cmd.remove('hydrogens')
    pymol.cmd.select('Pocket', f'byres {delete_extension_from_filename(ligand_filename)} around {distance}')
    pymol.cmd.save(output_pocket_path_full, 'Pocket')
    pymol.cmd.delete('all')


def generate_complex(ligand_path_full, pocket_path_full, output_complex_path_full):
    """
    Generates complex from the given ligand and computed pocket and saves it to a file.

    Args:
        ligand_path_full - full path to the input ligand file (including its name), must be a .pdb file!
        pocket_path_full - full path to the input pocket file (including its name), must be a .pdb file!
        output_complex_path_full - full path to the output complex file (including its name)
    """


    # # if the input ligand is not in a .pdb format - convert it to .pdb and change the path to the ligand
    # ligand_filename = extract_filename_from_full_path(ligand_path_full)
    # input_ligand_format = extract_format_from_filename(ligand_filename)
    # if input_ligand_format != 'pdb':
    #     # extract format of the ligand file
    #     ligand_foler = os.path.dirname(ligand_path_full)
    #     new_ligand_path_full = os.path.join(ligand_foler, delete_extension_from_filename(ligand_filename) + ".pdb")
    #     # os.system(f'obabel {ligand_input_path} -O {ligand_path} -d')
    #     obConversion = openbabel.OBConversion()
    #     obConversion.SetInAndOutFormats(input_ligand_format, "pdb")
    #     mol = openbabel.OBMol()
    #     obConversion.ReadFile(mol, ligand_path_full)
    #     obConversion.WriteFile(mol, new_ligand_path_full)
    #     ligand_path_full = new_ligand_path_full # change the path to the ligan after conversion

    assert ligand_path_full.endswith(".pdb"), "Ligand file must be a .pdb file!"
    assert pocket_path_full.endswith(".pdb"), "Pocket file must be a .pdb file!"

    ligand = Chem.MolFromPDBFile(ligand_path_full, sanitize=False, removeHs=False)
    if ligand == None:
        raise RuntimeError(f"Unable to process ligand: {ligand_path_full}")

    pocket = Chem.MolFromPDBFile(pocket_path_full, sanitize=False, removeHs=False)
    if pocket == None:
        raise RuntimeError(f"Unable to process pocket:  {ligand_path_full}")

    complex = (ligand, pocket)
    with open(output_complex_path_full, 'wb') as f:
        pickle.dump(complex, f)

def main(
    complex_name,
    db_path,
    base_correlation_name="_gnina_docked_cc.txt",
    threshold_correlation=0.6,
    base_protein_name="_protein_processed.pdb",
    distance=5, 
    is_main_log=True,
    main_log_filename="log.txt",
    main_log_path=os.path.join(os.getcwd(), "complex_main_logs"),

):
    """
    The main function that generates complexes (docked ligand + pocket) for those docking poses that passed
    cross-correlation check. Writes full paths to the complexes that passed cross-correlation check in a separate
    file: db_path/{complex_name}/{complex_name}_corr{threshold_correlation}_passed_complexes.txt for convenient future
    processing.

    Args:
        complex_name - name (id) of the protein-ligand complex
        db_path - path to the database with the complexes' data
        base_correlation_name - base name for the file with docking correlations (used to construct the full name)
        threshold_correlation - threshold value for the cross-correlation check
        base_protein_name - base name for the protein file (used to construct the full name)
        distance - threshold distance to generate the pocket
        is_main_log - whether to write logs for the main function
        main_log_filename - name of the log file for the main function
        main_log_path - path to the folder where main log files will be stored
    """

    try:
        complex_folder = os.path.join(db_path, complex_name) # folder with the complex data

        # list to store full paths to the docked ligands that passed cross-correlation check
        passed_docked_ligand_paths = []

        # read correlations and filter out docked ligands
        correlations_path_full = os.path.join(complex_folder, complex_name + base_correlation_name)
        with open(correlations_path_full, "r") as corr_file:
            for line in corr_file:
                ligand_path, corr = line.strip().split()
                corr = float(corr)
                if corr >= threshold_correlation:
                    passed_docked_ligand_paths.append(ligand_path)
        if len(passed_docked_ligand_paths) == 0:
            raise RuntimeError(f"No docked ligands passed cross-correlation check!")
        
        # find protein file in the database
        protein_path_full = find_pdb_protein_file_in_db(
            complex_name, db_path, base_protein_name=base_protein_name
        )

        # generate pockets and complexes for the filtered ligands and write paths to generated complexes
        # to a file (for convenient future processing)
        passed_complexes_path_file = os.path.join(
            complex_folder,
            f"{complex_name}_corr{threshold_correlation}_passed_complexes.txt"
        )
        with open(passed_complexes_path_file, "w") as passed_complexes_file:
            for passed_lignad_path in passed_docked_ligand_paths:
                lignad_fname = extract_filename_from_full_path(passed_lignad_path)
                output_pocket_path_full = os.path.join(
                    complex_folder,
                    f"{complex_name}_pocket_distance{distance}_{lignad_fname}"
                )
                generate_pocket(
                    passed_lignad_path,
                    protein_path_full,
                    output_pocket_path_full,
                    distance=distance
                )
                output_complex_path_full = os.path.join(
                    complex_folder,
                    f"{complex_name}_complex_distance{distance}_{delete_extension_from_filename(lignad_fname)}.rdkit"
                )
                generate_complex(
                    passed_lignad_path,
                    output_pocket_path_full,
                    output_complex_path_full
                )
                passed_complexes_file.write(f"{output_complex_path_full}\n")

        if is_main_log:
            log(
                f"Successfully generated complexes for {complex_name}. Check them in: {passed_complexes_path_file}.",
                status="INFO",
                log_filename=main_log_filename,
                log_path=main_log_path,
            )

    except Exception as e:
        if is_main_log:
            log(
                f"Failed to generate complexes for {complex_name}: {e}",
                status="ERROR",
                log_filename=main_log_filename,
                log_path=main_log_path,
            )

if __name__ == '__main__':

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

    # read script's arguments
    opts, args = getopt.getopt(sys.argv[1:], "s:e:p:")
    start_index = None  # start index for the complex names
    end_index = None  # end index for the complex names
    n_proc = None  # number of processes for multiprocessing (Pool.starmap)
    for opt, arg in opts:
        if opt == "-s":
            start_index = int(arg)
        elif opt == "-e":
            end_index = int(arg)
        elif opt == "-p":
            n_proc = int(arg)

    assert start_index is not None, "Start index is None after arugment's parsing"
    assert end_index is not None, "End index is None after arugment's parsing"
    assert n_proc is not None, "Number of processes is None after arugment's parsing"

    # apply start and end indice to the molecule names
    complex_names = complex_names[start_index : (end_index + 1)]
    n_complexes = len(complex_names)
    print(f"Computing good maps for {n_complexes} complexes....")

    # create log folders
    # NOTE: If is_..._log == True, it's important to create the folders cause the code won't work otherwise!
    is_main_log = True
    main_log_path = None
    if is_main_log:
        main_log_filename = (
            datetime.now(timezone.utc).strftime("%d_%m_%Y_%H.%M.%S.%f")[:-2]
            + "_main_log.txt"
        )
        main_log_path = os.path.join(os.getcwd(), "complexes_main_logs")
        create_folder(main_log_path)

    # specify keyword arguments for the main function
    base_correlation_name = "_gnina_docked_cc.txt"
    threshold_correlation = 0.6
    base_protein_name = "_protein_processed.pdb"
    distance = 5
    main_kwargs = {
        "base_correlation_name": base_correlation_name,
        "threshold_correlation": threshold_correlation,
        "base_protein_name": base_protein_name,
        "distance": distance, 
        "is_main_log": is_main_log,
        "main_log_filename": main_log_filename,
        "main_log_path": main_log_path,
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

    with Pool(n_proc) as pool:
        result = pool.starmap(apply_args_and_kwargs, starmap_args_iter)
        for r in result:
            pass