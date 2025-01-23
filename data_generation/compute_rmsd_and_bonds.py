import numpy as np
import os
import sys
import getopt
from itertools import repeat
from multiprocessing.pool import Pool
from rdkit.Chem import Lipinski

# append repo path to sys for convenient imports
sys.path.append(os.path.dirname(os.path.dirname(os.path.realpath(__file__))))

from utils import (
    read_molecule,
    find_pdb_ligand_file_in_db,
    delete_extension_from_filename,
    read_complexes_names_from_file,
    apply_args_and_kwargs,
)

def compute_rmsd(mol1, mol2):
    """
    Computes RMSD between two RDKit molecules.

    Args:
        mol1, mol2 - two RDKit molecules

    Returns:
        rmsd - computed value of the rmsd
    """

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
        output_path_full, 
        base_ligand_filename="_lignad.pdb",
        base_conformers_filename="_ligand_delatoms_conformers_nconfs3_generation_mode_gnina_docking.pdb",
    ):
    """
    The main function to compute rmsd between generated conformers and reference ligand structure. Writes computed
    rmsd, as well as, the number of rotatable bonds to an output file.

    Args:
        complex_name - name (id) of the protein-ligand complex
        db_path - path to the database with the complexes' data
        output_path_full - full path to the output file where rmsd and number of bonds will be written (same file for all complexes)
        base_ligand_name - base name for the ligand file (used to construct the full name)
        base_conformers_filename - base name for the file with conformers (used to construct the full name)
    """
    try:
        # find ligand file in the database and read it as an RDKit molecule
        ligand_path_full = find_pdb_ligand_file_in_db(complex_name, db_path, base_ligand_name=base_ligand_filename)
        ligand_mol = read_molecule(ligand_path_full, remove_Hs=True, sanitize=True)

        # read input conformers file and split it to sepearate files for each conformer
        conformers_path_full = find_pdb_ligand_file_in_db(complex_name, db_path, base_ligand_name=base_conformers_filename)
        with open(conformers_path_full, "r") as conf_file:
            lines = []
            temp_files_created = []
            count = 0
            for line in conf_file:
                if line == "\n":
                    continue
                lines.append(line)
                if line.strip() == "END":
                    count += 1
                    temp_file_path = os.path.join(db_path, complex_name, f"{complex_name}{delete_extension_from_filename(base_conformers_filename)}_temp_{count}.pdb")
                    with open(temp_file_path, "w") as temp_file:
                        temp_file.writelines(lines)
                    temp_files_created.append(temp_file_path)
                    lines = []
        
        # convert splitted conformers to RDKit molecules
        conformers_list = []
        for temp_file_path in temp_files_created:
            mol = read_molecule(temp_file_path, remove_Hs=True, sanitize=True)
            conformers_list.append(mol)
        
        # compute rmsd and number of rotatable bonds
        rmsd_list = []
        for conformer_mol in conformers_list:
            rmsd = compute_rmsd(ligand_mol, conformer_mol)
            rmsd_list.append(rmsd)
        avg_rmsd = sum(rmsd_list) / len(rmsd_list)
        rot_bonds = Lipinski.NumRotatableBonds(ligand_mol)

        # write rmsd and number of bonds to the output file
        with open(output_path_full, "a") as out_file:
            out_file.write(f"{complex_name} {avg_rmsd} {rot_bonds}\n")
        print(f"Computed rmsd for {complex_name}")
        
        # clear temporary files (splitted conformers' files)
        for temp_file_path in temp_files_created:
            os.remove(temp_file_path)
            
    except Exception as e:
        print(f"Failed to compute rmsd for complex {complex_name}: {e}")


    
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
        db_path + os.path.sep + "PDB_IDs_with_rdkit_length_less_than_24A.csv"
    )
    complex_names = read_complexes_names_from_file(complex_names_csv)
    n_complexes = len(complex_names)

    output_path_full = os.path.join(os.getcwd(), "rmsd_bonds.txt") # path to the output file where rmsd and number of bonds will be written

    base_ligand_filename="_ligand.pdb"
    base_conformers_filename = "_ligand_delatoms_conformers_nconfs3_generation_mode_gnina_docking.pdb"

    main_kwargs = {
        "base_ligand_filename": base_ligand_filename,
        "base_conformers_filename": base_conformers_filename,
    } # keyword arguments for the main() function

    # run computations in several Processes to speed up
    # generate args and kwargs iterables for Pool.starmap()
    main_args_iter = zip(
        complex_names,
        repeat(db_path, n_complexes),
        repeat(output_path_full, n_complexes)
    )  # iterable for positional arguments
    main_kwargs_iter = repeat(main_kwargs, n_complexes) # iterable for keyword arguments
    starmap_args_iter = zip(
        repeat(main, n_complexes), main_args_iter, main_kwargs_iter
    )  # combined positional and keyword arguments for Pool.starmap()

    with Pool(num_process) as pool:
        pool.starmap(apply_args_and_kwargs, starmap_args_iter)