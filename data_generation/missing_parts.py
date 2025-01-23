import os
import sys
import numpy as np
from subprocess import PIPE

# append repo path to sys for convenient imports
sys.path.append(os.path.dirname(os.path.dirname(os.path.realpath(__file__))))

from utils import (
    extract_filename_from_full_path,
)


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

