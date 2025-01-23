import os
import sys
import numpy as np
from subprocess import PIPE
from rdkit import Chem
from rdkit.Chem import AllChem

# append repo path to sys for convenient imports
sys.path.append(os.path.dirname(os.path.dirname(os.path.realpath(__file__))))

from utils import (
    read_molecule,
    save_all_conformers_to_pdb,
    delete_extension_from_filename,
    extract_filename_from_full_path,
    group_conformers_to_single_file,
    create_folder,
    split_sdf_file_to_pdbs,
)
from data_generation.Docking import gnina_docking

def rdkit_conformers(
    ligand_path_full,
    n_confs,
    conformers_path=os.path.join(os.getcwd(), "raw_conformers_data"),
    use_small_ring_torsions=False,
    prune_rms_tresh=1.0,
    random_seed=0xF00D,
):
    """
    Generates conformers using RDKit library and saves them to separate pdb files.

    Args:
        ligand_path_full - full path to the ligand file (icnluding its name), usually a .pdb file
        n_confs - number of conformers to generate
        conformers_path - path tp the folder where generated conformers will be written
        use_small_ring_torsions - whether to include additional small ring torsion potentials
        prune_rms_tresh - threshold value for RMSD pruning of the conformers
        random_seed - seed for conformers generation algorithm (for reproducibility)

    Returns:
        n_confs_generated - actual number of generated conformers (may be different from n_confs)
        conformer_path_full_list - list of full paths to the generated conformers
    """

    # read the ligand file
    mol = read_molecule(ligand_path_full, remove_Hs=False, sanitize=False)
    original_mol = Chem.Mol(mol)  # copy original molecule for future allignment

    # specify parameters for conformers generation
    params = AllChem.ETKDGv3()
    params.useSmallRingTorsions = use_small_ring_torsions
    params.pruneRMsThresh = prune_rms_tresh
    params.randomSeed = random_seed
    params.useRandomCoords = False

    # generate conformers
    AllChem.EmbedMultipleConfs(mol, numConfs=n_confs, params=params)

    # Align conformers to the original molecule
    conformer_ids = [x.GetId() for x in mol.GetConformers()]
    for conf_id in conformer_ids:
        AllChem.AlignMol(mol, original_mol, prbCid=conf_id, refCid=0)

    # write conformers to separate pdb files
    ligand_fname = extract_filename_from_full_path(
        ligand_path_full
    )  # extract name of the input ligand file
    base_conformer_filename = delete_extension_from_filename(ligand_fname) + ".pdb"
    n_confs_generated, conformer_path_full_list = save_all_conformers_to_pdb(
        mol, base_conformer_filename, pdb_path=conformers_path
    )

    return n_confs_generated, conformer_path_full_list


def gnina_conformers(
    ligand_path_full,
    protein_path_full,
    ligand_density_path_full,
    n_confs,
    conformers_path=os.path.join(os.getcwd(), "raw_conformers_data"),
    box_extension=5.0,
    path_to_gnina=os.path.join(os.getcwd(), "gnina"),
):
    """
    Generates conformers using GNINA docking and saves them to separate pdb files.

    Args:
        ligand_path_full - full path to the ligand file (icnluding its name), usually a .pdb file
        protein_path_full - full path to the protein file (icnluding its name), usually a .pdb file
        ligand_density_path_full - full path to the ligand density map which is used for the docking
        box size calculation
        n_confs - number of conformers to generate
        conformers_path - path to the folder where generated conformers will be written
        box_extension - used to compute the box size for docking, the box size is ligand_size + box_extension
        path_to_gnina - path to the gnina executable

    Returns:
        n_confs_generated - actual number of generated conformers (may be different from n_confs)
        conformer_path_full_list - list of full paths to the generated conformers

    """

    # generate output .sdf file for gnina docking
    ligand_fname = extract_filename_from_full_path(
        ligand_path_full
    )  # extract name of the input ligand file
    sdf_conformers_filename = (
        f"docked_confs_{delete_extension_from_filename(ligand_fname)}.sdf"
    )
    sdf_conformers_path_full = os.path.join(conformers_path, sdf_conformers_filename)

    # generate conformers with docking
    p = gnina_docking(
        ligand_path_full,
        protein_path_full,
        ligand_density_path_full,
        sdf_conformers_path_full,
        n_confs,
        box_extension=box_extension,
        path_to_gnina=path_to_gnina,
    )

    _, stderr = p.communicate()

    if p.returncode != 0 or (stderr and "open babel warning" not in stderr.decode().lower()):
        raise RuntimeError(f"Failed to perform docking: {stderr.decode().lower()}")

    # write conformers to separate pdb files
    base_conformer_filename = delete_extension_from_filename(ligand_fname) + ".pdb"
    n_confs_generated, conformer_path_full_list = split_sdf_file_to_pdbs(
        sdf_conformers_path_full,
        base_conformer_filename,
        pdb_path=conformers_path,
        remove_Hs=False,
    )

    return n_confs_generated, conformer_path_full_list


def generate_conformers(
    ligand_path_full,
    protein_path_full,
    ligand_density_path_full,
    n_confs,
    output_confs_path_full,
    generation_mode="gnina_docking",
    conformers_kwargs={},
    temp_data=os.path.join(os.getcwd(), "temp"),
):
    """
    Generates conformers, filters them out by cross-correlation threshold and saves
    filtered conformers into a single file (for convenient density map generation).

    Args:
        ligand_path_full - full path to the ligand file (icnluding its name), usually a .pdb file
        protein_path_full - full path to the protein file (icnluding its name), usually a .pdb file
        ligand_density_path_full - full path to the ligand density map which is used for the docking
        box size calculation
        n_confs - number of conformers to generate
        output_confs_path_full - full path to the output file with all conformers (icnluding its name), usually a .pdb file
        generation_mode - specifies mode for conformers generation (e.g. with docking or RDKit)
        conformers_kwargs - keyword parameters for conformers generation
        temp_data - folder where all temporary data for this function will be stored
        density_resolution - desired resolution of the map (in Angstrom)
        n_box - number of points for the map's cubic box
        grid_spacing - grid spacing for the cubic map
        is_chimeraX_log - should we write logs for ChimeraX scripts
        chimeraX_log_path - path to the folder where ChimeraX's log file will be stored
        (excluding the file's name which will be created automatically)
        chimeraX_script_path - path to the folder with the python scripts for ChimeraX
        threshold_correlation - threshold correlation value for cross-correlation check
        not_found_corr_value - the value we use if the correlation for a particular conformer wasn't found/computed
        chimeraX_output_base_filename - base name of the file to store output from ChimeraX script
        (needed to read correlation value for a particular pose)
        chimeraX_output_path - path to the folder where output ChimeraX files will be stored
        clear_chimeraX_output - whether to clear ChimeraX output files
        write_corrs_to_file - whether to write computed correlations to a file
        corrs_path_full - full path to the file (including its name) where computed correlations will be written
        (required only if write_corrs_to_file == True)

    Returns:
        n_confs_generated - number of generated conformers (NOTE: may be different from n_confs)
    """

    # create a folder for conformers generation
    conformers_path = os.path.join(temp_data, "raw_conformers_data")
    create_folder(conformers_path)

    # generate conformers
    generation_mode = generation_mode.strip().lower()
    match generation_mode:
        case "gnina_docking":
            n_confs_generated, conformer_path_full_list = gnina_conformers(
                ligand_path_full,
                protein_path_full,
                ligand_density_path_full,
                n_confs,
                conformers_path=conformers_path,
                **conformers_kwargs,
            )

        case "rdkit":
            n_confs_generated, conformer_path_full_list = rdkit_conformers(
                ligand_path_full,
                n_confs,
                conformers_path=conformers_path,
                **conformers_kwargs,
            )

        case default:
            raise RuntimeError(
                f"Conformers generation mode {generation_mode} is not implemented!"
            )
        
    if n_confs_generated == 0:
        raise RuntimeError(
            f"No conformers were generated using {generation_mode}!"
        )
    
    # group generated conformers into a single file for convenient processing
    group_conformers_to_single_file(
        conformer_path_full_list, output_confs_path_full, delete_input=False
    )

    return n_confs_generated

if __name__ == "__main__":

    ligand_path_full = os.path.join(
        os.getcwd(), "maps_generated_with_model", "2zv9_ligand.pdb"
    )

    protein_path_full = os.path.join(
        os.getcwd(), "maps_generated_with_model", "2zv9_protein_processed.pdb"
    )

    ligand_density_path_full = os.path.join(
        os.getcwd(), "maps_generated_with_model", "2zv9_ligand_res_2.0_gridpsace_0.5_nbox_48_size_24A_label.mrc"
    )

    lignad_fname = extract_filename_from_full_path(
        ligand_path_full
    )  # extract ligand's file name

    # create a folder to store temporary data
    temp_data = os.path.join(
        os.getcwd(), f"temp_{delete_extension_from_filename(lignad_fname)}"
    )
    create_folder(temp_data)

    # specify parameters for conformers generation
    n_confs = 3
    output_confs_path_full = os.path.join(
        os.getcwd(), "maps_generated_with_model", f"{delete_extension_from_filename(lignad_fname)}_conformers.pdb"
    )
    generation_mode = "gnina_docking"
    box_extension = 1.0
    conformers_kwargs = {
        "box_extension": box_extension,
        "path_to_gnina": os.path.join(os.getcwd(), "gnina"),
    }

    # generate conformers
    generate_conformers(
        ligand_path_full,
        protein_path_full,
        ligand_density_path_full,
        n_confs,
        output_confs_path_full,
        generation_mode=generation_mode,
        conformers_kwargs=conformers_kwargs,
        temp_data=temp_data,
    )