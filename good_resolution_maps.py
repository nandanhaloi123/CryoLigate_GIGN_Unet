import os
import getopt
import sys
from itertools import repeat
from datetime import datetime, timezone
from multiprocessing.pool import Pool
from utils import (
    find_pdb_ligand_file_in_db,
    compute_density_map_in_chimeraX,
    read_complexes_names_from_file,
    create_folder,
    apply_args_and_kwargs,
    log,
)


def main(
    complex_name,
    db_path,
    base_ligand_name="_ligand.pdb",
    is_main_log=True,
    main_log_filename="log.txt",
    main_log_path=os.path.join(os.getcwd(), "docking_main_logs"),
    density_resolution=1.0,
    n_box=16,
    is_chimeraX_log=True,
    chimeraX_log_path=os.path.join(os.getcwd(), "chimeraX_logs"),
    chimeraX_script_path=os.path.join(os.getcwd(), "chimeraX_scripts"),
):
    """
    The main function for computing good resolution density map.

    Args:
        complex_name - name (id) of the protein-ligand complex
        db_path - path tp the database with the complexes' data
        base_ligand_name - base name for the ligand file (used to construct the full name)
        is_main_log - whether to write logs for the main function
        main_log_filename - name of the log file for the main function
        main_log_path - path to the folder where main log files will be stored
        density_resolution - desired resolution of the map (in Angstrom) that we generate for docking poses
        n_box - number of points for the map's cubic box
        is_chimeraX_log - should we write logs for ChimeraX scripts
        chimeraX_log_path - path to the folder where ChimeraX's log file will be stored
        (excluding the file's name which will be created automatically)
        chimeraX_script_path - path to the folder with the python scripts for ChimeraX
    """
    try:
        # find ligand file in the database
        ligand_path_full = find_pdb_ligand_file_in_db(
            complex_name, db_path, base_ligand_name=base_ligand_name
        )

        # compute the map
        complex_folder = os.path.join(db_path, complex_name) # folder with the complex data
        output_density_path_full = os.path.join(
            complex_folder,
            f"{complex_name}_ligand_res{density_resolution}_boxed_{n_box}A.mrc"
        )
        p = compute_density_map_in_chimeraX(
            ligand_path_full,
            output_density_path_full,
            density_resolution=density_resolution,
            n_box=n_box,
            is_log=is_chimeraX_log,
            log_path=chimeraX_log_path,
            script_path=chimeraX_script_path
        )
        _, stderr = p.communicate() # catch ChimeraX errors if any
        if p.returncode != 0 or stderr:
            raise RuntimeError(f"Failed to compute density map: {stderr}")

        if is_main_log:
            log(
                f"Successfully computed good resolution map for {complex_name}. Check the map in: {output_density_path_full}.",
                status="INFO",
                log_filename=main_log_filename,
                log_path=main_log_path,
            )

    except Exception as e:
        if is_main_log:
            log(
                f"Failed to compute good resolution map for {complex_name}: {e}",
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
        main_log_path = os.path.join(os.getcwd(), "good_maps_main_logs")
        create_folder(main_log_path)

    # specify keyword arguments for the main function
    base_ligand_name = "_ligand.pdb"
    density_resolution = 1.0
    n_box = 16
    chimeraX_script_path = os.path.join(os.getcwd(), "chimeraX_scripts")
    main_kwargs = {
        "base_ligand_name": base_ligand_name,
        "is_main_log": is_main_log,
        "main_log_filename": main_log_filename,
        "main_log_path": main_log_path,
        "density_resolution": density_resolution,
        "n_box": n_box,
        "is_chimeraX_log": is_chimeraX_log,
        "chimeraX_log_path": chimeraX_log_path,
        "chimeraX_script_path": chimeraX_script_path,
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