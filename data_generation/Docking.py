

######### Perform docking ############### 
# from vina import Vina
import os
import sys
import getopt
import shutil
import pandas as pd
import mrcfile
import numpy as np
from itertools import repeat
from multiprocessing.pool import Pool
from datetime import datetime, timezone
from subprocess import Popen, PIPE, DEVNULL
from tqdm import tqdm
from rdkit import Chem
from rdkit.Chem import AllChem
from pymol import cmd

# append repo path to sys for convenient imports
sys.path.append(os.path.dirname(os.path.dirname(os.path.realpath(__file__))))

from utils import (
    split_sdf_file_to_pdbs, 
    create_folder, 
    compute_mol_map_correlation_in_chimeraX, 
    extract_correlation_from_chimera_file,
    read_density_data_mrc,
    calculate_center_of_mass_of_density,
    calculate_ligand_size,
    find_pdb_ligand_file_in_db,
    find_pdb_protein_file_in_db,
    find_mrc_density_file_in_db,
    log,
    apply_args_and_kwargs,
    read_complexes_names_from_file,

)


def gnina_docking(    
    ligand_path_full,
    protein_path_full,
    ligand_density_path_full,
    output_path_full,
    n_pos,
    box_extension=5.0,
    path_to_gnina=os.getcwd() + os.path.sep + "gnina" 
):
    """
    Performs docking using gnina software. Runs gnina commands through a subprocess.
    
    Args:
        ligand_path_full - full path to the ligand file (icnluding its name), must be a .pdb file!
        protein_path_full - full path to the protein file (icnluding its name), usually a .pdb file
        ligand_density_path_full - full path to the file with ligand's density (icnluding its name), 
        must be an .mrc file!
        output_path_full - full path to the output file with generated poses (icnluding its name), 
        usually a .pdb file
        n_pos - number of poses to generate in docking
        box_extension - used to compute the box size for docking, the box size is ligand_size + box_extension
        path_to_gnina - path to the gnina executable

    Returns:
        p - object of the Popen class corresponding to the gnina's subprocess
    """

    assert ligand_density_path_full.endswith(".mrc"), "Error in gnina docking: Density file must be an .mrc file!"
    assert ligand_path_full.endswith(".pdb"), "Error in gnina docking: Ligand file must be a .pdb file!"

    density_map, header, voxel_size = read_density_data_mrc(ligand_density_path_full) # TODO: replace with low resolution maps 
    origin = header.origin
    origin = np.array(origin.tolist())
    
    center_of_mass = calculate_center_of_mass_of_density(density_map, voxel_size=voxel_size, origin=origin) 

    center_of_mass_x = center_of_mass[0]
    center_of_mass_y = center_of_mass[1]
    center_of_mass_z = center_of_mass[2]

    # Calculate the molecular size
    mol = Chem.MolFromPDBFile(ligand_path_full, sanitize=False, removeHs=False)
    ligand_size = calculate_ligand_size(mol)

    # compute box size
    box_size = ligand_size + box_extension

    p = Popen(
        [
            path_to_gnina,
            "-r", protein_path_full,
            "-l", ligand_path_full,
            "--center_x", str(center_of_mass_x),
            "--center_y", str(center_of_mass_y),
            "--center_z", str(center_of_mass_z),
            "--size_x", str(box_size),
            "--size_y", str(box_size),
            "--size_z", str(box_size),
            "--num_modes", str(n_pos),
            "-o", output_path_full,
            "--cpu", str(1)
        ],
        stdout=DEVNULL,
        stderr=PIPE,
    )

    return p


def filter_docked_poses_by_correlation(
    docked_path_full_list,
    target_density_path_full,
    threshold_correlation=0.6,
    not_found_corr_value=0.0,
    chimeraX_output_base_filename="output.txt",
    chimeraX_output_path=os.getcwd() + os.path.sep + "chimeraX_output",
    density_resolution=1.0,
    n_box=16,
    grid_spacing=0.5,
    is_log=False,
    log_path=os.getcwd() + os.path.sep + "chimeraX_logs",
    script_path=os.getcwd() + os.path.sep + "chimeraX_scripts",
    clear_chimeraX_output=True,
    write_corrs_to_file=True,
    corrs_path_full=os.getcwd() + os.path.sep + "docking_correlations.txt"
):
    """
    Filters out generated docking poses by given threshold correlation.
    The correlation is computed between the density map of each provided pose and
    the target density that was used in docking. For a particular pose, if the correlation 
    is less than the threshod value then this pose is considered to be inappropriate. 

    Args:
        docked_path_full_list - list of full paths to the input docked poses (including their filenames)
        target_density_path_full - full path to the target denisty file (including its filename)
        threshold_correlation - threshold correlation value
        not_found_corr_value - the value we use if the correlation for a particular pose wasn't found/computed
        chimeraX_output_base_filename - base name of the file to store output from ChimeraX script 
        (needed to read correlation value for a particular pose)
        chimeraX_output_path - path to the folder where output ChimeraX files will be stored
        density_resolution - desired resolution of the maps we generate for the molecules (in Angstrom)
        n_box - number of points for the cubic box
        grid_spacing - grid spacing for the cubic map
        is_log - should we write logs for ChimeraX scripts
        log_path - path to the folder where the log file will be stored (excluding the file's name which will be created automatically)
        script_path - path to the folder with the python script for ChimeraX (excluding its name)
        clear_chimeraX_output - whether to clear ChimeraX output files
        write_corrs_to_file - whether to write computed correlations to a file
        corrs_path_full - full path to the file (including its name) where computed correlations will be written
        (required only if write_corrs_to_file == True)

    Returns:
        corrs - list with computed correlations for the docked poses
        n_appr - number of the appropriate poses (based on the threshold correlation)
        appr_dock_path_full_list - list of full paths to the appropriate docked poses
    """

    # create folder for ChimeraX output files if doesn't exist
    create_folder(chimeraX_output_path)

    corrs = [] # list to store computed correlations
    n_appr = 0 # number of appropriate docked poses (those that passed the correlation threshold)
    appr_dock_path_full_list = [] # list to store full paths to the appropriate docked poses

    # compute correlations for the provided poses and filter them by threshold value
    for i in range(len(docked_path_full_list)):
        docked_path_full = docked_path_full_list[i]

        # construct full path to the output ChimeraX file
        output_filename =  f"corr_pos_{i+1}_{chimeraX_output_base_filename}"
        output_path_full = chimeraX_output_path + os.path.sep + output_filename 

        # compute correlation inside ChimeraX
        corr = not_found_corr_value
        with open(output_path_full, "w") as output_file: # ChimeraX output with correlation will be written here
            p = compute_mol_map_correlation_in_chimeraX(
                docked_path_full,
                target_density_path_full,
                output_file,
                density_resolution=density_resolution,
                n_box=n_box,
                grid_spacing=grid_spacing,
                is_log=is_log,
                log_path=log_path,
                script_path=script_path,
            )
            _, stderr = p.communicate() # catch errors from ChimeraX subprocess

        if (p.returncode == 0) and (not stderr):
            try:
                corr = extract_correlation_from_chimera_file(
                    output_path_full, 
                    not_found_value=not_found_corr_value
                )
            except Exception as e:
                print(f"Failed to extract correlation from Chimera output file: {e}")

        # store computed correlation
        corrs.append(corr)

        # filter poses by the threshold correlation 
        if corr >= threshold_correlation:
            n_appr += 1
            appr_dock_path_full_list.append(docked_path_full)
    
    # write computed correlation to file if needed
    if write_corrs_to_file:
        with open(corrs_path_full, "w") as file:
            for i in range(len(docked_path_full_list)):
                docked_path_full = docked_path_full_list[i]
                corr = corrs[i]
                file.write(f"{docked_path_full} {corr}\n") 
    
    # clear ChimeraX output if needed
    if clear_chimeraX_output:
        shutil.rmtree(chimeraX_output_path)

    return corrs, n_appr, appr_dock_path_full_list


def main(
    complex_name,
    db_path,
    base_ligand_name="_ligand.pdb",
    base_protein_name="_protein_processed.pdb",
    base_density_name="_ligand.mrc",
    n_pos=10,
    box_extension=5.0,
    path_to_gnina=os.getcwd() + os.path.sep + "gnina", 
    is_main_log=True,
    main_log_filename="log.txt",
    main_log_path=os.path.join(os.getcwd(), "docking_main_logs"),
    write_corrs_to_file=True,
    threshold_correlation=0.6,
    not_found_corr_value=0.0,
    density_resolution=1.0,
    n_box=32,
    grid_spacing=0.5,
    is_chimeraX_log=True,
    chimeraX_log_path=os.path.join(os.getcwd(), "chimeraX_logs"),
    chimeraX_script_path=os.path.join(os.getcwd(), "chimeraX_scripts"),
    clear_chimeraX_output=True,
):
    """
    The main function that generates docked poses. If required, computes cross correlations for them and 
    writes the correlations to a file.

    Args:
        complex_name - name (id) of the protein-ligand complex
        db_path - path tp the database with the complexes' data
        base_ligand_name - base name for the ligand file (used to construct the full name)
        base_protein_name - base name for the protein file (used to construct the full name)
        base_density_name - base name for the density file (used to construct the full name)
        n_pos - number of poses to generate in docking
        box_extension - used to compute the box size for docking, the box size is ligand_size + box_extension
        path_to_gnina - path to the gnina executable
        is_main_log - whether to write logs for the main function
        main_log_filename - name of the log file for the main function
        main_log_path - path to the folder where main log files will be stored
        write_corrs_to_file - whether to write computed correlations to a file 
        (will be written in db_path/{complex_name}/{complex_name}_gnina_docked_cc.txt)
        threshold_correlation - threshold correlation value for cross-correlation check
        not_found_corr_value - the value we use if the correlation for a particular pose wasn't found/computed
        density_resolution - desired resolution of the map (in Angstrom) that we generate for docking poses
        n_box - number of points for the map's cubic box
        grid_spacing - grid spacing for the cubic map
        is_chimeraX_log - should we write logs for ChimeraX scripts
        chimeraX_log_path - path to the folder where ChimeraX's log file will be stored
        (excluding the file's name which will be created automatically)
        chimeraX_script_path - path to the folder with the python scripts for ChimeraX
        clear_chimeraX_output - whether to clear ChimeraX output files that are used for extracting cross-correlation values
    """
    try:
        # find ligand file in the database
        ligand_path_full = find_pdb_ligand_file_in_db(
            complex_name, db_path, base_ligand_name=base_ligand_name
        )

        # find protein file in the database
        protein_path_full = find_pdb_protein_file_in_db(
            complex_name, db_path, base_protein_name=base_protein_name
        )

        # find density file in the database
        density_path_full = find_mrc_density_file_in_db(
            complex_name, db_path, base_density_name=base_density_name
        )

        # do docking
        # generate output .sdf file for gnina docking
        complex_folder = os.path.join(db_path, complex_name)
        sdf_docking_filename = (
            f"{complex_name}_gnina_docked.sdf"
        )
        sdf_docking_path_full = os.path.join(complex_folder, sdf_docking_filename)

        # run gnina inside a subporcess
        p = gnina_docking(
            ligand_path_full,
            protein_path_full,
            density_path_full,
            sdf_docking_path_full,
            n_pos,
            box_extension=box_extension,
            path_to_gnina=path_to_gnina,
        )
        _, stderr = p.communicate() # catch gnina errors if any
        if p.returncode != 0 or (stderr and "open babel warning" not in stderr.decode().lower()):
            raise RuntimeError(f"{stderr}")

        # split output sdf file to separate pdb files (for convenient cross-correlation check)
        base_docking_filename = "gnina_docked.pdb"
        n_poses_saved, poses_path_full_list = split_sdf_file_to_pdbs(
            sdf_docking_path_full,
            base_docking_filename,
            pdb_path=complex_folder,
            remove_Hs=False,
        )
        if n_poses_saved == 0:
            raise RuntimeError(f"No docking poses were saved for {complex_name}!")
        
        # compute cross-correlations and write them to a file if required
        if write_corrs_to_file:
            chimeraX_output_base_filename = "output.tmp"
            chimeraX_output_path = os.path.join(os.getcwd(), f"{complex_name}_chimeraX_output")
            create_folder(chimeraX_output_path) # create a folder for ChimeraX output if it doesn't exist
            corrs_path_full = os.path.join(
                complex_folder,
                f"{complex_name}_gnina_docked_cc.txt"
            ) # full path to the fil where correlations will be written
            filter_docked_poses_by_correlation(
                poses_path_full_list,
                density_path_full,
                threshold_correlation=threshold_correlation,
                not_found_corr_value=not_found_corr_value,
                chimeraX_output_base_filename=chimeraX_output_base_filename,
                chimeraX_output_path=chimeraX_output_path,
                density_resolution=density_resolution,
                n_box=n_box,
                grid_spacing=grid_spacing,
                is_log=is_chimeraX_log,
                log_path=chimeraX_log_path,
                script_path=chimeraX_script_path,
                clear_chimeraX_output=clear_chimeraX_output,
                write_corrs_to_file=write_corrs_to_file,
                corrs_path_full=corrs_path_full,
            )

        if is_main_log:
            log(
                f"Successfully finished docking for {complex_name}. Check the data in: {complex_folder}.",
                status="INFO",
                log_filename=main_log_filename,
                log_path=main_log_path,
            )

    except Exception as e:
        if is_main_log:
            log(
                f"Failed to compute docking for {complex_name}: {e}",
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

    # load molecule names
    complex_names_csv = (
        db_path + os.path.sep + "PDB_IDs_with_rdkit_length_less_than_24A_nbonds_less_than_10.csv"
    )
    complex_names = read_complexes_names_from_file(complex_names_csv)

    # apply start and end indices to the molecule names
    complex_names = complex_names[start_index : (end_index + 1)]
    n_complexes = len(complex_names)
    print(f"Computing docking for {n_complexes} complexes....")

    # create log folders
    is_chimeraX_log = False
    chimeraX_log_path = None
    if is_chimeraX_log:
        chimeraX_log_path = os.path.join(os.getcwd(), "chimeraX_logs")
        create_folder(chimeraX_log_path)
    is_main_log = True
    main_log_path = None
    if is_main_log:
        main_log_filename = (
            datetime.now(timezone.utc).strftime("%d_%m_%Y_%H.%M.%S.%f")[:-2]
            + "_main_log.txt"
        )
        main_log_path = os.path.join(os.getcwd(), "docking_main_logs")
        create_folder(main_log_path)

    # specify keyword arguments for the main function
    base_ligand_name = "_ligand.pdb"
    base_protein_name = "_protein_processed.pdb"
    base_density_name = "_nconfs3_genmode_gnina_docking_boxextens1.0_res4.0_delprob0.0_low_resolution_forward_model.mrc"
    n_pos = 3
    box_extension = 1.0
    path_to_gnina = os.getcwd() + os.path.sep + "gnina" 
    write_corrs_to_file = True
    threshold_correlation = 0.6
    not_found_corr_value = 0.0
    density_resolution = 4.0
    n_box = 48
    grid_spacing = 0.5
    chimeraX_script_path = os.path.join(os.getcwd(), "chimeraX_scripts")
    clear_chimeraX_output = True
    main_kwargs = {
        "base_ligand_name": base_ligand_name,
        "base_protein_name": base_protein_name,
        "base_density_name": base_density_name,
        "n_pos": n_pos,
        "box_extension": box_extension,
        "path_to_gnina": path_to_gnina, 
        "is_main_log": is_main_log,
        "main_log_filename": main_log_filename,
        "main_log_path": main_log_path,
        "write_corrs_to_file": write_corrs_to_file,
        "threshold_correlation": threshold_correlation,
        "not_found_corr_value": not_found_corr_value,
        "density_resolution": density_resolution,
        "n_box": n_box,
        "grid_spacing": grid_spacing,
        "is_chimeraX_log": is_chimeraX_log,
        "chimeraX_log_path": chimeraX_log_path,
        "chimeraX_script_path": chimeraX_script_path,
        "clear_chimeraX_output": clear_chimeraX_output,
    }

    # run computations in several Processes to speed up
    # generate args and kwargs iterables for Pool.starmap()
    main_args_iter = zip(
        complex_names,
        repeat(db_path, n_complexes),
    )  # iterable for positional arguments
    main_kwargs_iter = repeat(main_kwargs, n_complexes) # iterable for keyword arguments
    starmap_args_iter = zip(
        repeat(main, n_complexes), main_args_iter, main_kwargs_iter
    )  # combined positional and keyword arguments for Pool.starmap()

    with Pool(n_proc) as pool:
        result = pool.starmap(apply_args_and_kwargs, starmap_args_iter)
        for r in result:
            pass

    


