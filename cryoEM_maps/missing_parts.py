import os
import sys
import numpy as np
from subprocess import PIPE

# append repo path to sys for convenient imports
sys.path.append(os.path.dirname(os.path.dirname(os.path.realpath(__file__))))

from utils import (
    compute_density_map_in_chimera,
    read_density_data_mrc,
    delete_extension_from_filename,
    extract_filename_from_full_path,
    create_folder,
)

# NOTE: we can try this version of the function in the future
# def random_delete_atoms_from_pdb_file(
#     pdb_path_full,
#     delatoms_molecule_path=os.getcwd() + os.path.sep + "delatoms_molecule_data",
#     delete_prob=0.2,
# ):
#     """
#     Randomly deletes atoms from the given pdb file to simulate the case when
#     some parts of the molecule are completely missed.

#     Params:
#     pdb_path_full - full path to the input .pdb file (icnluding its name)
#     delatoms_molecule_path - path to the directory where files with some deleted atoms are
#     stored
#     delete_prob - probability for deleting an atom

#     Returns:
#     delatoms_pdb_path_full - full path to the new .pdb file with some atoms deleted
#     """

#     assert delete_prob >= 0.0, "Probability of atoms deleting should be >= 0"

#     # exctract pdf filename from the input full path
#     pdb_filename = extract_filename_from_full_path(pdb_path_full)

#     assert pdb_filename.endswith(".pdb"), "pdb filename must end with .pdb!"

#     # create directory for the molecule files with some deleted atoms if it doesn't exist
#     create_folder(delatoms_molecule_path)

#     # read input lines
#     with open(pdb_path_full, "r") as input_file:
#         input_lines = input_file.readlines()

#     # split input file lines into those related and unrelated to atoms
#     atoms_lines = [] # list to store indices of lines related to atoms
#     non_atoms_lines = set() # set to store indices of lines unrelated to atoms
#     for i, input_line in enumerate(input_lines):
#         if input_line.startswith("ATOM") or input_line.startswith("HETATM"):
#             atoms_lines.append(i)
#         else:
#             non_atoms_lines.add(i)

#     # randomly delete lines related to atoms
#     n_atom_lines_to_delete = int(len(atoms_lines) * delete_prob) # the number of atom lines to delete
#     # generate indices of atoms' lines that should be kept
#     atom_lines_to_keep = set(np.random.choice(atoms_lines, size=len(atoms_lines) - n_atom_lines_to_delete, replace=False))

#     # write a new file with some atoms deleted
#     delatoms_pdb_path_full = delatoms_molecule_path + os.path.sep + "del_atoms_" + pdb_filename # full path to the pdb file with some deleted atoms

#     with open(delatoms_pdb_path_full, "w") as output_file:
#         for i, input_line in enumerate(input_lines):
#             if i in non_atoms_lines or i in atom_lines_to_keep:
#                 output_file.write(input_line)

#     return delatoms_pdb_path_full


def random_delete_atoms_from_pdb_file(
    input_pdb_path_full,
    delatoms_pdb_path_full,
    delete_prob=0.2,
):
    """
    Randomly deletes atoms from the given pdb file to simulate the case when
    some parts of the molecule are completely missed.

    Args:
        input_pdb_path_full - full path to the input .pdb file (icnluding its name)
        delatoms_pdb_path_full - full path to the output file with some deleted atoms (icnluding its name)
        delete_prob - probability for deleting an atom
    """

    # exctract pdb filename from the input full path
    input_pdb_filename = extract_filename_from_full_path(input_pdb_path_full)
    assert input_pdb_filename.endswith(".pdb"), "pdb filename must end with .pdb!"

    # exctract pdb filename from the output full path
    delatoms_pdb_filename = extract_filename_from_full_path(delatoms_pdb_path_full)
    assert delatoms_pdb_filename.endswith(".pdb"), "pdb filename must end with .pdb!"

    with open(input_pdb_path_full, "r") as input_file:
        with open(delatoms_pdb_path_full, "w") as output_file:
            for input_line in input_file:
                if input_line.startswith("ATOM") or input_line.startswith("HETATM"):
                    # if the current line corresponds to an atom, decide whether this atom should
                    # be deleted by using uniform random number between 0 and 1.
                    random_val = np.random.rand()
                    if random_val > delete_prob:
                        output_file.write(input_line)
                else:
                    output_file.write(input_line)


if __name__ == "__main__":

    # construct input paths
    input_filename = "1c3b_ligand.pdb"
    input_path = os.getcwd() + os.path.sep + "raw_molecule_data"
    input_path_full = input_path + os.path.sep + input_filename

    # delete some atoms from the input pdb file and write a new file
    delete_prob = 0.2
    del_pdb_path = random_delete_atoms_from_pdb_file(
        input_path_full, delete_prob=delete_prob
    )  # returns full path to new pdb file

    # run python script inside Chimera to compute density map
    density_resolution = 3.5  # resolution of the density map (in Angstrom)
    density_path = (
        os.getcwd() + os.path.sep + "density_maps"
    )  # path to the folder where density will be stored
    density_filename = delete_extension_from_filename(input_filename) + ".mrc"
    create_folder(
        density_path
    )  # create a folder for the density files if it doesn't exist
    density_path_full = density_path + os.path.sep + density_filename

    # create a folder for Chimera logs if it doesn't exist
    chimera_log_path = os.getcwd() + os.path.sep + "chimera_logs"
    create_folder(chimera_log_path)

    p = compute_density_map_in_chimera(
        del_pdb_path,
        density_path_full,
        is_log=True,
        log_path=chimera_log_path,
        density_resolution=density_resolution,
        stderr_file=PIPE,
    )
    _, stderr = p.communicate()  # output from the subprocess

    # if the subprocess finished with errors, check logs in chimera_logs folder
    if p.returncode != 0 or stderr:
        print(stderr)
        print("Failed to compute density map. Check logs in chimera_logs folder.")

    # if the subprocess finished correctly, check obtained denisty map
    # in the densisty_maps folder
    else:
        density = read_density_data_mrc(density_path_full)
        print(density.shape)
        print(np.nonzero(density))
