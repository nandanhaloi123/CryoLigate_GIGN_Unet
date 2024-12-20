import os
import getopt
import sys
import numpy as np
import shutil
from itertools import repeat
from multiprocessing.pool import Pool
from datetime import datetime, timezone

# append repo path to sys for convenient imports
sys.path.append(os.path.dirname(os.path.dirname(os.path.realpath(__file__))))

from utils import (
    compute_density_map_in_chimeraX,
    log,
    create_folder,
    extract_filename_from_full_path,
    delete_extension_from_filename,
    find_pdb_ligand_file_in_db,
    find_pdb_protein_file_in_db,
    read_complexes_names_from_file,
    apply_args_and_kwargs,
)

from cryoEM_maps.internal_flexibility import (
    generate_conformers,
)

from cryoEM_maps.missing_parts import random_delete_atoms_from_pdb_file

def generate_low_resolution_density(
    ligand_path_full,
    protein_path_full,
    output_density_path_full,
    temp_data=os.path.join(os.getcwd(), "temp"),
    n_confs=10,
    generation_mode="docking",
    conformers_kwargs={},
    density_resolution=3.5,
    n_box=16,
    is_chimeraX_log=True,
    chimeraX_log_path=os.path.join(os.getcwd(), "chimeraX_logs"),
    chimeraX_script_path=os.path.join(os.getcwd(), "chimeraX_scripts"),
    threshold_correlation=0.6,
    not_found_corr_value=0.0,
    write_corrs_to_file=True,
    delete_prob=0.2,
):
    """
    Creates low resolution density map for the given ligand.

    Args:
        ligand_path_full - full path to the ligand file (icnluding its name), usually a .pdb file
        protein_path_full - full path to the protein file (icnluding its name), usually a .pdb file
        output_density_path_full - full path to the output density file (icnluding its name), usually a .mrc file
        temp_data - folder where all temporary data for this function will be stored
        n_confs - number of conformers to generate
        generation_mode - specifies mode for conformers generation (e.g. with docking or RDKit)
        conformers_kwargs - keyword parameters for conformers generation
        density_resolution - desired resolution of the map (in Angstrom)
        n_box - number of points for the map's cubic box
        is_chimeraX_log - should we write logs for ChimeraX scripts
        chimeraX_log_path - path to the folder where ChimeraX's log file will be stored
        (excluding the file's name which will be created automatically)
        chimeraX_script_path - path to the folder with the python scripts for ChimeraX
        threshold_correlation - threshold correlation value for cross-correlation check
        not_found_corr_value - the value we use if the correlation for a particular conformer wasn't found/computed
        write_corrs_to_file - whether to write computed correlations to a file (will be written in temp_data/docking_correlations.txt)
        delete_prob - probability for deleting an atom
    """
    # parameters for conformers generation
    lignad_fname = extract_filename_from_full_path(
        ligand_path_full
    )  # extract ligand's file name
    output_confs_path_full = os.path.join(
        temp_data, f"{delete_extension_from_filename(lignad_fname)}_conformers.pdb"
    )  # construct path to the output file with conformers
    chimeraX_output_base_filename = "output.tmp"
    chimeraX_output_path = os.path.join(temp_data, "chimeraX_output")
    create_folder(chimeraX_output_path)
    clear_chimeraX_output = False
    corrs_path_full = os.path.join(temp_data, "docking_correlations.txt")

    # generate conformers
    generate_conformers(
        ligand_path_full,
        protein_path_full,
        n_confs,
        output_confs_path_full,
        generation_mode=generation_mode,
        conformers_kwargs=conformers_kwargs,
        temp_data=temp_data,
        density_resolution=density_resolution,
        n_box=n_box,
        is_chimeraX_log=is_chimeraX_log,
        chimeraX_log_path=chimeraX_log_path,
        chimeraX_script_path=chimeraX_script_path,
        threshold_correlation=threshold_correlation,
        not_found_corr_value=not_found_corr_value,
        chimeraX_output_base_filename=chimeraX_output_base_filename,
        chimeraX_output_path=chimeraX_output_path,
        clear_chimeraX_output=clear_chimeraX_output,
        write_corrs_to_file=write_corrs_to_file,
        corrs_path_full=corrs_path_full,
    )

    # randomly delete atoms from conformers
    del_atoms_conformers_path_full = os.path.join(
        temp_data, f"del_atoms_confs_{delete_extension_from_filename(lignad_fname)}.pdb"
    )
    random_delete_atoms_from_pdb_file(
        output_confs_path_full,
        del_atoms_conformers_path_full,
        delete_prob=delete_prob,
    )

    # generate final density map for the conformers
    p = compute_density_map_in_chimeraX(
        del_atoms_conformers_path_full,
        output_density_path_full,
        density_resolution=density_resolution,
        n_box=n_box,
        is_log=is_chimeraX_log,
        log_path=chimeraX_log_path,
        script_path=chimeraX_script_path,
    )

    _, stderr = p.communicate()

    if p.returncode != 0 or stderr:
        raise RuntimeError(f"Failed to compute final density map: {stderr}")


def main(
    complex_name,
    db_path,
    output_density_path_full,
    base_ligand_name="_ligand.pdb",
    base_protein_name="_protein_processed.pdb",
    n_confs=10,
    generation_mode="docking",
    conformers_kwargs={},
    density_resolution=3.5,
    n_box=16,
    is_chimeraX_log=True,
    chimeraX_log_path=os.path.join(os.getcwd(), "chimeraX_logs"),
    chimeraX_script_path=os.path.join(os.getcwd(), "chimeraX_scripts"),
    threshold_correlation=0.6,
    not_found_corr_value=0.0,
    write_corrs_to_file=True,
    delete_prob=0.2,
    is_main_log=True,
    main_log_filename="log.txt",
    main_log_path=os.path.join(os.getcwd(), "maps_main_logs"),
    clear_temp=True,
):
    """
    The main function for generating low resolution density maps for the complexes in the database.

    Args:
        complex_name - name (id) of the protein-ligand complex
        db_path - path tp the database with the complexes' data
        output_density_path_full - full path to the output density file (icnluding its name), usually a .mrc file
        base_ligand_name - base name for the ligand file (used to construct the full name)
        base_protein_name - base name for the protein file (used to construct the full name)
        n_confs - number of conformers to generate
        generation_mode - specifies mode for conformers generation (e.g. with docking or RDKit)
        conformers_kwargs - keyword parameters for conformers generation
        density_resolution - desired resolution of the map (in Angstrom)
        n_box - number of points for the map's cubic box
        is_chimeraX_log - should we write logs for ChimeraX scripts
        chimeraX_log_path - path to the folder where ChimeraX's log file will be stored
        (excluding the file's name which will be created automatically)
        chimeraX_script_path - path to the folder with the python scripts for ChimeraX
        threshold_correlation - threshold correlation value for cross-correlation check
        not_found_corr_value - the value we use if the correlation for a particular conformer wasn't found/computed
        write_corrs_to_file - whether to write computed correlations to a file (will be written in temp_data/docking_correlations.txt)
        delete_prob - probability for deleting an atom
        is_main_log - whether to write logs for the main function
        main_log_filename - name of the log file for the main function
        main_log_path - path to the folder where main log files will be stored
        clear_temp - whether to delete temporary data when the function is done
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

        # create a folder for temporary data
        temp_data = os.path.join(os.getcwd(), f"{complex_name}_temp")
        create_folder(temp_data)

        # generate map
        generate_low_resolution_density(
            ligand_path_full,
            protein_path_full,
            output_density_path_full,
            temp_data=temp_data,
            n_confs=n_confs,
            generation_mode=generation_mode,
            conformers_kwargs=conformers_kwargs,
            density_resolution=density_resolution,
            n_box=n_box,
            is_chimeraX_log=is_chimeraX_log,
            chimeraX_log_path=chimeraX_log_path,
            chimeraX_script_path=chimeraX_script_path,
            threshold_correlation=threshold_correlation,
            not_found_corr_value=not_found_corr_value,
            write_corrs_to_file=write_corrs_to_file,
            delete_prob=delete_prob,
        )

        if is_main_log:
            log(
                f"Successfully computed density map for {complex_name}. Check the map in: {output_density_path_full}.",
                status="INFO",
                log_filename=main_log_filename,
                log_path=main_log_path,
            )

    except Exception as e:
        if is_main_log:
            log(
                f"Failed to compute density map for {complex_name}: {e}",
                status="ERROR",
                log_filename=main_log_filename,
                log_path=main_log_path,
            )

    finally:
        if clear_temp:  # delete the folder with temporary data if needed
            shutil.rmtree(temp_data)


if __name__ == "__main__":

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
    print(f"Computing low resolution density maps for {n_complexes} complexes....")

    # create log folders
    # NOTE: If is_..._log == True, it's important to create the folders cause the code won't work otherwise!
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
        main_log_path = os.path.join(os.getcwd(), "main_logs")
        create_folder(main_log_path)

    # specify keyword arguments for the main function
    base_ligand_name = "_ligand.pdb"
    base_protein_name = "_protein_processed.pdb"
    n_confs = 10
    generation_mode = "gnina_docking"
    conformers_kwargs = {
        "box_extension": 5.0,
        "path_to_gnina": os.path.join(os.getcwd(), "gnina"),
    }
    density_resolution = 3.5
    n_box = 16
    chimeraX_script_path = os.path.join(os.getcwd(), "chimeraX_scripts")
    threshold_correlation = 0.3
    not_found_corr_value = 0.0
    write_corrs_to_file = False
    delete_prob = 0.2
    clear_temp = True
    main_kwargs = {
        "base_ligand_name": base_ligand_name,
        "base_protein_name": base_protein_name,
        "n_confs": n_confs,
        "generation_mode": generation_mode,
        "conformers_kwargs": conformers_kwargs,
        "density_resolution": density_resolution,
        "n_box": n_box,
        "is_chimeraX_log": is_chimeraX_log,
        "chimeraX_log_path": chimeraX_log_path,
        "chimeraX_script_path": chimeraX_script_path,
        "threshold_correlation": threshold_correlation,
        "not_found_corr_value": not_found_corr_value,
        "write_corrs_to_file": write_corrs_to_file,
        "delete_prob": delete_prob,
        "is_main_log": is_main_log,
        "main_log_filename": main_log_filename,
        "main_log_path": main_log_path,
        "clear_temp": clear_temp,
    }

    # run computations in several Processes to speed up
    # generate args and kwargs iterables for Pool.starmap()
    outup_density_path_list = (
        os.path.join(
            db_path,
            complex_name,
            f"{complex_name}_nconfs{n_confs}_genmode_{generation_mode}_res{density_resolution}_nbox{n_box}_threshcorr{threshold_correlation}_delprob{delete_prob}_low_resolution_forward_model.mrc",
        )
        for complex_name in complex_names
    )
    main_args_iter = zip(
        complex_names,
        repeat(db_path, n_complexes),
        outup_density_path_list,
    )  # iterable for positional arguments
    main_kwargs_iter = repeat(main_kwargs, n_complexes)
    starmap_args_iter = zip(
        repeat(main, n_complexes), main_args_iter, main_kwargs_iter
    )  # combined positional and keyword arguments for Pool.starmap()

    with Pool(n_proc) as pool:
        result = pool.starmap(apply_args_and_kwargs, starmap_args_iter)
        for r in result:
            pass
