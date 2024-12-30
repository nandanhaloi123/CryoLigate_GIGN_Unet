import os
import getopt
import sys
from itertools import repeat
from multiprocessing.pool import Pool
from utils import delete_extension_from_filename, read_complexes_names_from_file, apply_args_and_kwargs


def main(complex_name, db_path):

    complex_folder = os.path.join(db_path, complex_name)
    txt_path = os.path.join(complex_folder, f"{complex_name}_complexes_map_less_noisy_bad_to_map2.0.txt")

    with open(txt_path, "w") as file:
        file.write(os.path.join(complex_folder, f"{complex_name}_5A.rdkit"))

    print(f"Created txt file: {txt_path}")
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
    n_complexes = len(complex_names)

    opts, args = getopt.getopt(sys.argv[1:], "p:")
    n_proc = None  # number of processes for multiprocessing (Pool.starmap)
    for opt, arg in opts:
        if opt == "-p":
            n_proc = int(arg)

    assert n_proc is not None, "Number of processes is None after arugment's parsing"

    main_args_iter = zip(
        complex_names,
        repeat(db_path, n_complexes),
    )  # iterable for positional arguments
    main_kwargs = {}
    main_kwargs_iter = repeat(main_kwargs, n_complexes)
    starmap_args_iter = zip(
        repeat(main, n_complexes), main_args_iter, main_kwargs_iter
    )  # combined positional and keyword arguments for Pool.starmap()

    with Pool(n_proc) as pool:
        result = pool.starmap(apply_args_and_kwargs, starmap_args_iter)
        for r in result:
            pass