import os
import numpy as np
import pickle
import torch
import mrcfile
from datetime import datetime, timezone
from rdkit import Chem
from subprocess import Popen, PIPE, DEVNULL


def normalize(x):
    return (x - x.min()) / (x.max() - x.min())

def create_dir(dir_list):
    assert  isinstance(dir_list, list) == True
    for d in dir_list:
        if not os.path.exists(d):
            os.makedirs(d)

def save_model_dict(model, model_dir, msg):
    model_path = os.path.join(model_dir, msg + '.pt')
    torch.save(model.state_dict(), model_path)
    print("model has been saved to %s." % (model_path))

def load_model_dict(model, ckpt):
    model.load_state_dict(torch.load(ckpt))

def del_file(path):
    for i in os.listdir(path):
        path_file = os.path.join(path,i)  
        if os.path.isfile(path_file):
            os.remove(path_file)
        else:
            del_file(path_file)

def write_pickle(filename, obj):
    with open(filename, 'wb') as f:
        pickle.dump(obj, f)

def read_pickle(filename):
    with open(filename, 'rb') as f:
        obj = pickle.load(f)
    return obj

class BestMeter(object):
    """Computes and stores the best value"""

    def __init__(self, best_type):
        self.best_type = best_type  
        self.count = 0      
        self.reset()

    def reset(self):
        if self.best_type == 'min':
            self.best = float('inf')
        else:
            self.best = -float('inf')

    def update(self, best):
        self.best = best
        self.count = 0

    def get_best(self):
        return self.best

    def counter(self):
        self.count += 1
        return self.count


class AverageMeter(object):
    """Computes and stores the average and current value"""

    def __init__(self):
        self.reset()

    def reset(self):
        self.val = 0
        self.avg = 0
        self.sum = 0
        self.count = 0

    def update(self, val, n=1):
        self.val = val
        self.sum += val * n
        self.count += n

    def get_average(self):
        self.avg = self.sum / (self.count + 1e-12)

        return self.avg

def run_script_inside_chimera(
    script_name,
    script_path=os.getcwd() + os.path.sep + "chimera_scripts",
    option_names=(),
    option_values=(),
    output_file=DEVNULL
):
    """
    Runs given python script inside Chimera through a subprocess (simulates running
    script from the terminal command:
    chimera --nogui --nostatus --script "name_of_script.py").
    NOTE: if some option in option_names doesn't have a value, an empty string ""
    must be provided to the option_values.
    Returns an object of the corresponding subprocess.

    Args:
        script_name - name of the python script
        script_path - path to the python script (excluding its name)
        option_names - names of the options provided to the script
        option_values - values for the options provided to the script
        output_file - file object where Chimera output will be written

    Returns:
        p - object of the Popen class corresponding to the Chimera's subprocess
    """

    assert script_name.endswith(
        ".py"
    ), "The name of the python script must end with .py"

    assert len(option_names) == len(
        option_values
    ), "Mismatch between number of options and number of values provided"

    # full path to the python script (including its name)
    script_path_full = script_path + os.path.sep + script_name

    # construct the line of the terminal command that describes the script
    script_line = script_path_full
    for i in range(len(option_names)):
        if option_values[i]:  # if option value is not an empty string
            option = f" {option_names[i]}" + f" {option_values[i]}"
        else:
            option = f" {option_names[i]}"
        script_line += option
    
    # construct full Chimera's command
    command = [
            "chimera",
            "--nogui",
            "--nostatus",
            "--script", script_line,
            ]   

    p = Popen(
        command,
        stdout=output_file,
        stderr=PIPE,
    )

    return p


def run_script_inside_chimeraX(
    script_name,
    script_path=os.getcwd() + os.path.sep + "chimeraX_scripts",
    option_names=(),
    option_values=(),
    output_file=DEVNULL,
):
    """
    Runs given python script inside ChimeraX through a subprocess (simulates running
    script from the terminal command:
    ChimeraX --nogui --nostatus --script "name_of_script.py").
    NOTE: if some option in option_names doesn't have a value, an empty string ""
    must be provided to the option_values.
    Returns an object of the corresponding subprocess.

    Args:
        script_name - name of the python script
        script_path - path to the python script (excluding its name)
        option_names - names of the options provided to the script
        option_values - values for the options provided to the script
        output_file - file object where ChimeraX output will be written

    Returns:
        p - object of the Popen class corresponding to the ChimeraX's subprocess
    """

    assert script_name.endswith(
        ".py"
    ), "The name of the python script must end with .py"

    assert len(option_names) == len(
        option_values
    ), "Mismatch between number of options and number of values provided"

    # full path to the python script (including its name)
    script_path_full = script_path + os.path.sep + script_name

    # construct the line of the terminal command that describes the script
    script_line = script_path_full
    for i in range(len(option_names)):
        if option_values[i]:  # if option value is not an empty string
            option = f" {option_names[i]}" + f" {option_values[i]}"
        else:
            option = f" {option_names[i]}"
        script_line += option
    
    # construct full Chimera's command
    command = [
            "ChimeraX",
            "--nogui",
            "--nostatus",
            "--script", script_line,
            ]   

    p = Popen(
        command,
        stdout=output_file,
        stderr=PIPE,
    )

    return p


def compute_density_map_in_chimera(
    molecule_path_full,
    density_path_full,
    density_resolution=1.0,
    is_log=False,
    log_path=os.getcwd() + os.path.sep + "chimera_logs",
    script_path=os.getcwd() + os.path.sep + "chimera_scripts",
):
    """
    Computes density map for the given molecule/conformer using Chimera.
    To achieve this, runs python script "chimera_density_map.py" inside Chimera through a
    subprocess.

    Args:
        molecule_path_full - full path to the input molecule file (including its name)
        density_path_full - full path to the output density file (including its name)
        density_resolution - desired resolution of the map (in Angstrom)
        is_log - should we write logs for Chimera scripts
        log_path - path to the folder where the log file will be stored (excluding the file's name which will be created automatically)
        script_path - path to the folder with the python script for Chimera (excluding its name)

    Returns:
        p - object of the Popen class corresponding to the Chimera's subprocess
    """

    # name of the python script to run inside Chimera
    scipt_name = "chimera_density_map.py"

    # options for the python script
    option_names = ["-i", "-r", "-o"]
    option_values = [molecule_path_full, density_resolution, density_path_full]
    if is_log:
        option_names.append("-l")
        option_values.append(log_path)

    # returns object corresponding to the subprocess
    p = run_script_inside_chimera(
        scipt_name,
        script_path=script_path,
        option_names=option_names,
        option_values=option_values,
    )

    return p


def compute_density_map_in_chimeraX(
    molecule_path_full,
    density_path_full,
    density_resolution=1.0,
    n_box=16,
    grid_spacing=0.5,
    is_log=False,
    log_path=os.getcwd() + os.path.sep + "chimeraX_logs",
    script_path=os.getcwd() + os.path.sep + "chimeraX_scripts",
):
    """
    Computes density map for the given molecule/conformer using ChimeraX. 
    The computed map has a cubic box with the specified number of points (n_box).
    To achieve this, runs python script "chimeraX_density_map.py" inside ChimeraX through a
    subprocess.

    Args:
        molecule_path_full - full path to the input molecule file (including its name)
        density_path_full - full path to the output density file (including its name)
        density_resolution - desired resolution of the map (in Angstrom)
        n_box - number of points for the cubic box
        grid_spacing - grid spacing for the cubic map
        is_log - should we write logs for ChimeraX scripts
        log_path - path to the folder where the log file will be stored (excluding the file's name which will be created automatically)
        script_path - path to the folder with the python script for ChimeraX (excluding its name)

    Returns:
        p - object of the Popen class corresponding to the ChimeraX's subprocess
    """

    # name of the python script to run inside ChimeraX
    scipt_name = "chimeraX_density_map.py"

    # options for the python script
    option_names = ["-i", "-r", "-b", "-g", "-o"]
    option_values = [molecule_path_full, density_resolution, n_box, grid_spacing, density_path_full]
    if is_log:
        option_names.append("-l")
        option_values.append(log_path)

    # returns object corresponding to the subprocess
    p = run_script_inside_chimeraX(
        scipt_name,
        script_path=script_path,
        option_names=option_names,
        option_values=option_values,
    )

    return p


def compute_density_map_on_grid_in_chimeraX(
    molecule_path_full,
    target_map_path_full,
    output_density_path_full,
    density_resolution=1.0,
    is_log=False,
    log_path=os.getcwd() + os.path.sep + "chimeraX_logs",
    script_path=os.getcwd() + os.path.sep + "chimeraX_scripts",
):
    """
    Computes density map for the given molecule/conformer on the grid of the target map using ChimeraX. 
    To achieve this, runs python script "chimeraX_density_map_on_grid.py" inside ChimeraX through a
    subprocess.

    Args:
        molecule_path_full - full path to the input molecule file (including its name)
        target_map_path_full - full path to the target map which grid will be used to compute 
        the molecule's map
        output_density_path_full - full path to the molecule's output density file (including its name)
        density_resolution - desired resolution of the map (in Angstrom)
        is_log - should we write logs for ChimeraX scripts
        log_path - path to the folder where the log file will be stored (excluding the file's name which will be created automatically)
        script_path - path to the folder with the python script for ChimeraX (excluding its name)

    Returns:
        p - object of the Popen class corresponding to the ChimeraX's subprocess
    """

    # name of the python script to run inside ChimeraX
    scipt_name = "chimeraX_density_map_on_grid.py"

    # options for the python script
    option_names = ["-i", "-t", "-r", "-o"]
    option_values = [molecule_path_full, target_map_path_full, density_resolution, output_density_path_full]
    if is_log:
        option_names.append("-l")
        option_values.append(log_path)

    # returns object corresponding to the subprocess
    p = run_script_inside_chimeraX(
        scipt_name,
        script_path=script_path,
        option_names=option_names,
        option_values=option_values,
    )

    return p


def compute_mol_map_correlation_in_chimera(
    molecule_path_full,
    target_density_path_full,
    output_file,
    density_resolution=1.0,
    is_log=False,
    log_path=os.getcwd() + os.path.sep + "chimera_logs",
    script_path=os.getcwd() + os.path.sep + "chimera_scripts",
):
    """
    Computes correlation between density map of a given molecule and target density map using Chimera.
    To achieve this, runs python script "chimera_mol_density_corr.py" inside Chimera through a
    subprocess.

    Args:
        molecule_path_full - full path to the input molecule file (including its name)
        target_density_path_full - full path to the target density file (including its name)
        output_file - file object where Chimera output will be written. NOTE: here it's necessarry to provide this
        file since the correlation will be written there
        density_resolution - desired resolution of the map we generate for the molecule (in Angstrom)
        is_log - whether we should write logs for Chimera scripts
        log_path - path to the folder where the log file will be stored (excluding the file's name which will be created automatically)
        script_path - path to the folder with the python script for Chimera (excluding its name)

    Returns:
        p - object of the Popen class corresponding to the Chimera's subprocess
    """

    # name of the python script to run inside Chimera
    scipt_name = "chimera_mol_density_corr.py"

    # options for the python script
    option_names = ["-i", "-r", "-t"]
    option_values = [molecule_path_full, density_resolution, target_density_path_full]
    if is_log:
        option_names.append("-l")
        option_values.append(log_path)

    # returns object corresponding to the subprocess
    p = run_script_inside_chimera(
        scipt_name,
        script_path=script_path,
        option_names=option_names,
        option_values=option_values,
        output_file=output_file
    )

    return p


def compute_mol_map_correlation_in_chimeraX(
    molecule_path_full,
    target_density_path_full,
    output_file,
    density_resolution=1.0,
    n_box=16,
    grid_spacing=0.5,
    is_log=False,
    log_path=os.getcwd() + os.path.sep + "chimeraX_logs",
    script_path=os.getcwd() + os.path.sep + "chimeraX_scripts",
):
    """
    Computes correlation between density map of a given molecule and target density map using ChimeraX.
    To achieve this, runs python script "chimeraX_mol_density_corr.py" inside ChimeraX through a
    subprocess. NOTE: Assumes that both density maps are cubic with the specified number of points (n_box).

    Args:
        molecule_path_full - full path to the input molecule file (including its name)
        target_density_path_full - full path to the target density file (including its name)
        output_file - file object where ChimeraX output will be written. NOTE: here it's necessarry to provide this
        file since the correlation will be written there
        density_resolution - desired resolution of the map we generate for the molecule (in Angstrom)
        n_box - number of points for the cubic box
        grid_spacing - grid spacing for the cubic map
        is_log - whether we should write logs for ChimeraX scripts
        log_path - path to the folder where the log file will be stored (excluding the file's name which will be created automatically)
        script_path - path to the folder with the python script for ChimeraX (excluding its name)

    Returns:
        p - object of the Popen class corresponding to the ChimeraX's subprocess
    """

    # name of the python script to run inside ChimeraX
    scipt_name = "chimeraX_mol_density_corr.py"

    # options for the python script
    option_names = ["-i", "-r", "-b", "-g", "-t"]
    option_values = [molecule_path_full, density_resolution, n_box, grid_spacing, target_density_path_full]
    if is_log:
        option_names.append("-l")
        option_values.append(log_path)

    # returns object corresponding to the subprocess
    p = run_script_inside_chimeraX(
        scipt_name,
        script_path=script_path,
        option_names=option_names,
        option_values=option_values,
        output_file=output_file
    )

    return p

def fit_into_map_in_chimera(
    input_path_full,
    target_density_path_full,
    output_path_full,
    output_chimera_file,
    is_log=False,
    log_path=os.getcwd() + os.path.sep + "chimera_logs",
    script_path=os.getcwd() + os.path.sep + "chimera_scripts",
):
    """
    Fits given input molecule/density into the target density using Chimera.
    To achieve this, runs python script "chimera_fit_into_map.py" inside Chimera through a
    subprocess.


    Args:
        input_path_full - full path to the input molecule/density file (including its name)
        target_density_path_full - full path to the target density file (including its name)
        output_path_full - full path to the output file (including its name) where fitted molecule/density will be written
        output_chimera_file - file object where Chimera output will be written. NOTE: here it's necessarry to provide this
        file since fitting scores will be written there
        is_log - whether we should write logs for Chimera scripts
        log_path - path to the folder where the log file will be stored (excluding the file's name which will be created automatically)
        script_path - path to the folder with the python script for Chimera (excluding its name)

    Returns:
        p - object of the Popen class corresponding to the Chimera's subprocess
    """

    # name of the python script to run inside Chimera
    scipt_name = "chimera_fit_into_map.py"

    # options for the python script
    option_names = ["-i", "-t", "-o"]
    option_values = [input_path_full, target_density_path_full, output_path_full]
    if is_log:
        option_names.append("-l")
        option_values.append(log_path)

    # returns object corresponding to the subprocess
    p = run_script_inside_chimera(
        scipt_name,
        script_path=script_path,
        option_names=option_names,
        option_values=option_values,
        output_file=output_chimera_file
    )

    return p

def extract_correlation_from_chimera_file(output_path_full, not_found_value=0.0):
    """
    Extracts correlation value from the given Chimera output file. 
    NOTE: works only for both Chimera and ChimeraX.

    Args:
        output_path_full - full path to Chimera output file (including its name)
        not_found_value - the value which will be returned if the correlation value
        is not found in the file

    Returns:
        correlation - correlation value   
    """

    # assign initial value, it will be returned if the correlation value
    # is not found in the file
    correlation = not_found_value

    with open(output_path_full, "r") as file:
        for line in file:
            line = line.strip().lower()
            if line.startswith("correlation"):
                start_ind = line.find("=") # index of the equl sign
                if start_ind == -1: # if equal sign wasn't found, go to the next line
                    continue
                start_ind += 1 # junp over the whitespace
                end_ind = start_ind # index of the comma separating correlation value
                while line[end_ind] != ",":
                    end_ind += 1
                correlation = float(line[start_ind: end_ind])
                break

    return correlation


# NOTE: for now, the function commented out below shouldn't be used in the code

# def prepare_raw_pdb_file(
#     raw_pdb_filename,
#     raw_molecule_path=os.getcwd() + os.path.sep + "raw_molecule_data",
#     molecule_path=os.getcwd() + os.path.sep + "molecule_data",
# ):
#     """
#     Prepares a raw .pdb molecule file to feed into Chimera: deletes all unnecessary rows.

#     Params:
#     raw_pdb_filename - name of the raw input .pdb file
#     raw_molecule_path - path to the directory where raw molecule files are stored
#     molecule_path - path to the directory where processed (prepared) molecule files are stored

#     Returns:
#     prepared_pdb_filename - name of the prepared_pdb_file
#     """

#     assert raw_pdb_filename.endswith(".pdb"), "pdb filename must end with .pdb!"

#     # create directory for prepared molecule files if it doesn't exist
#     create_folder(molecule_path)

#     # construct full path to the input raw pdb file
#     raw_pdb_path_full = raw_molecule_path + os.path.sep + raw_pdb_filename

#     # construct full path to the prepared pdb file
#     prepared_pdb_filename = "prepared_" + raw_pdb_filename
#     pdb_path_full = molecule_path + os.path.sep + prepared_pdb_filename

#     with open(raw_pdb_path_full, "r") as raw_file:
#         with open(pdb_path_full, "w") as prepared_file:
#             for raw_line in raw_file:
#                 # Keep only the lines corresponding to atoms
#                 if raw_line.startswith("ATOM") or raw_line.startswith("HETATM"):
#                     prepared_file.write(raw_line)

#     return prepared_pdb_filename


def read_density_data_txt(
    density_filename,
    density_path=os.getcwd() + os.path.sep + "density_maps",
):
    """
    Reads density map data from the given .txt file and converts it to numpy array.

    Args:
        density_filename - name of the file with density data
        density_path - path to the directory where density file is stored

    Returns:
        density - numpy array with density data
    """

    assert density_filename.endswith(".txt"), "txt filename must end with .txt!"

    # full path to the density file (including its name)
    path_full = density_path + os.path.sep + density_filename

    # fetch dimensions of the density matrix from the file
    dims = tuple()
    with open(path_full, "r") as file:
        first_line = file.readline()
        first_line = first_line.split(":")[-1]  # dimensions are written here
        first_line = first_line[1:-2]  # remove \n and parentheses
        dims = tuple(map(int, first_line.split(", ")))

    density = np.loadtxt(path_full)
    density = density.reshape(dims)

    return density

def delete_extension_from_filename(filename):
    """
    Deletes extension (symbols after .) from the given filename.

    Args:
        filename - name of the given file

    Returns:
        name of the file without extension
    """

    return ".".join(filename.split(".")[:-1])


def extract_filename_from_full_path(path_full):
    """
    Extracts filename from the full (including the name) path to a file.

    Args:
        full_path - full path to the file

    Returns:
        filename extracted from the full path
    """

    return path_full.split(os.path.sep)[-1]


def extract_folder_from_full_path(path_full):
    """
    Extracts path to the folder from the full (including the name) path to a file.

    Args:
        full_path - full path to the file

    Returns:
        path to the folder extracted from the full path
    """

    return os.path.sep.join(path_full.split(os.path.sep)[:-1])


def extract_format_from_filename(filename):
    """
    Gets file's format from the given filename.

    Args:
        filename - name of the file

    Returns:
        format of the file
    """

    return filename.split(".")[-1]


def read_molecules_from_sdf(
    sdf_path_full,
    n_mols=1,
    remove_Hs=False,
    sanitize=False,
):
    """
    Reads molecules data from the given .sdf file using RDKit library.

    Args:
        sdf_path_full - full path to the input .sdf file (including its name)
        n_mols - expected number of molecules to read from the file
        remove_Hs - whether to remove hydrogen atoms from the molecule when reading
        sanitize - whether to sanitize molecule using RDKit

    Returns:
        mols - list with RDKit molecule objects obtained from the given file
    """

    assert sdf_path_full.endswith(".sdf"), "sdf filename must end with .sdf!"

    mols = []  # list to store RDKit molecule objects

    with Chem.SDMolSupplier(sdf_path_full, sanitize=sanitize, removeHs=remove_Hs) as suppl:
        for mol in suppl:
            if mol is None:
                print(
                    f"Failed to read one of the molecules from the sdf file in: {sdf_path_full}"
                )
                continue
            mols.append(mol)

    assert (
        len(mols) == n_mols
    ), f"Expected {n_mols} molecules from the sdf file, but read {len(mols)} molecules. File path: {sdf_path_full}"

    return mols


def read_molecule_from_pdb(
    pdb_path_full,
    remove_Hs=False,
    sanitize=False,
):
    """
    Reads one molecule data from the given .pdb file using RDKit library.
    NOTE: even if you have several molecules in the given file, RDKit will
    read only the first one.

    Args:
        pdb_path_full - full path to the input .pdb file (including its name)
        remove_Hs - whether to remove hydrogen atoms from the molecule when reading
        sanitize - whether to sanitize molecule using RDKit
    Returns:
        mol - RDKit molecule object obtained from the given file
    """

    assert pdb_path_full.endswith(".pdb"), "pdb filename must end with .pdb!"

    mol = Chem.MolFromPDBFile(pdb_path_full, sanitize=sanitize, removeHs=remove_Hs)

    return mol


def read_molecule(
    path_full, remove_Hs=True, sanitize=False,
):
    """
    Reads one molecule data from the given file and converts it to an RDKit molecule object.

    Args:
        path_full - full path to the input molecule file
        remove_Hs - whether to remove hydrogen atoms from the molecule when reading
        sanitize - whether to sanitize molecule using RDKit

    Returns:
        mol - RDKit molecule object obtained from the given file
    """

    # exctract name of the file from full path
    filename = extract_filename_from_full_path(path_full)

    # extract format of the name
    file_format = extract_format_from_filename(filename)

    match file_format:
        case "sdf":
            mols = read_molecules_from_sdf(
                path_full, n_mols=1, remove_Hs=remove_Hs, sanitize=sanitize,
            )
            return mols[0]

        case "pdb":
            mol = read_molecule_from_pdb(
                path_full, remove_Hs=remove_Hs, sanitize=sanitize,
            )
            return mol

        case default:
            raise RuntimeError(
                f"Reading from .{file_format} molecule files is not implemented!"
            )


def save_one_conformer_to_pdb(
    mol,
    conf_id,
    pdb_filename,
    pdb_path=os.getcwd() + os.path.sep + "raw_conformers_data",
):
    """
    Saves specified (by conf_id) conformer of the given RDKit molecule to a .pdb file.

    Args:
        mol - RDKit molecule object
        conf_id - id of the molecule's conformer to save. -1 indicates the last conformer
        pdb_filename - name of the pdb file
        pdb_path - path to the folder where the pdb file should stored (excluding its name)

    Returns:
        pdb_path_full - full path to the conformer's pdb file (including its name)
    """

    assert pdb_filename.endswith(".pdb"), "pdb filename must end with .pdb!"

    # construct full path to the pdb file (including its name)
    pdb_path_full = pdb_path + os.path.sep + pdb_filename

    Chem.MolToPDBFile(
        mol,
        pdb_path_full,
        confId=conf_id,
    )

    return pdb_path_full


def save_all_conformers_to_pdb(
    mol,
    base_pdb_filename,
    pdb_path=os.getcwd() + os.path.sep + "raw_conformers_data",
):
    """
    Saves all conformers of the given RDKit molecule to .pdb files.
    Each conformer is saved in a separate file.

    Args:
        mol - RDKit molecule object
        base_pdb_filename - base name of the output pdb files. "conformer_{conformer_id}" will
        be added to this name
        pdb_path - path to the folder where the pdb files should be stored

    Returns:
        n_conf - number of saved conformers
        pdb_path_full_list - list with full paths to the conformers' pdb files
    """

    assert base_pdb_filename.endswith(".pdb"), "pdb filename must end with .pdb!"

    pdb_path_full_list = []  # list with full paths to the output conformers' files

    # iterate over all available conformers
    conf_ids = [x.GetId() for x in mol.GetConformers()]
    for conf_id in conf_ids:
        # construct the name of the pdb file for the given conformer
        pdb_filename = f"conf_{conf_id + 1}_" + base_pdb_filename
        pdb_path_full = save_one_conformer_to_pdb(
            mol, conf_id, pdb_filename, pdb_path=pdb_path
        )
        pdb_path_full_list.append(pdb_path_full)

    n_conf = len(pdb_path_full_list)

    return n_conf, pdb_path_full_list


def group_conformers_to_single_file(
    path_full_list,
    group_path_full,
    delete_input=True,
):
    """
    Groups multiple conformer files into a single file. The paths to the input conformers files are
    provided in the path_full_list.

    Args:
        path_full_list - list with full paths to the input conformers' files to group (including their names)
        group_path_full - full path to the output group file (including its name)
        delete_input - whether to delete input conformers' files
    """
    # extract group filename from full path
    group_filename = extract_filename_from_full_path(group_path_full)

    # exctract format of the group file
    group_file_format = extract_format_from_filename(group_filename)

    # check if the format of each input file match with the group file format
    for path_full in path_full_list:
        input_filename = extract_filename_from_full_path(path_full)
        input_file_format = extract_format_from_filename(input_filename)
        assert (
            input_file_format == group_file_format
        ), f"Format of the group file ({group_file_format}) isn't the same as input file's format ({input_file_format}). Change input file {path_full} or group file {group_path_full}."

    # group input files into a sinlge output file
    with open(group_path_full, "w") as group_file:
        for path_full in path_full_list:
            with open(path_full, "r") as input_file:
                group_file.writelines(input_file)

            group_file.write("\n")

            if delete_input:
                os.remove(path_full)


def rescale_density_map(
        density_path_full, 
        rescaled_path_full, 
        box_size=16,
    ):
    """
    Rescales given density map such that it fits the given box size. 
    Achieves this by running relion software commands through a subprocess.

    Args:
        density_path_full - full path to the input denisty map i.e. map to rescale (including its name)
        rescaled_path_full - full path to the output file with rescaled density
        box_size - size of the box (the box to which we rescale)

    Returns:
        p - object of the Popen class corresponding to the relion's subprocess
    """
  
    # run relion's command through a subpprocess
    p = Popen(
    [
        "relion_image_handler",
        "--i", density_path_full,
        "--new_box", str(box_size),
        "--o", rescaled_path_full,
    ],
    stderr=PIPE,
    )
    
    return p


def align_center_mass_of_densities(to_align_path_full, target_path_full):
    """
    Shifts given density map (to_align_path_full) such that its
    center of mass matches the one of target density map (target_path_full).
    Both files must be .mrc files!

    Args:
        to_align_path_full - full path to the density we align (including file's name)  
        target_path_full - full path to the target density (including file's name)   
    """
    
    assert to_align_path_full.endswith(".mrc"), "Input mrc filename must end with .mrc!"
    assert target_path_full.endswith(".mrc"), "Target mrc filename must end with .mrc!"

    # read data of the map to align
    density_map_to_align, header_to_align, voxel_size_to_align = read_density_data_mrc(to_align_path_full)
    origin_to_align = header_to_align.origin
    origin_to_align = np.array(origin_to_align.tolist())

    # for the map to align compute coordinates of center of mass
    center_of_mass_to_align = calculate_center_of_mass_of_density(
                                                                density_map_to_align, 
                                                                voxel_size=voxel_size_to_align, 
                                                                origin=origin_to_align
                                                                ) 
    
    # read data of the target map
    density_map_target, header_target, voxel_size_target = read_density_data_mrc(target_path_full)
    origin_target = header_target.origin
    origin_target = np.array(origin_target.tolist())

    # for the target map compute indices of center of mass (center of mass in voxel space)
    center_of_mass_target = calculate_center_of_mass_of_density(
                                                                density_map_target,
                                                                voxel_size=voxel_size_target, 
                                                                origin=origin_target
                                                                ) 
    
    # calculate the difference between centers of mass
    diff = center_of_mass_target - center_of_mass_to_align

    with mrcfile.open(to_align_path_full, 'r+') as mrc:
        mrc.header.origin.x += diff[0]
        mrc.header.origin.y += diff[1]
        mrc.header.origin.z += diff[2]

def align_density_maps(
        to_align_path_full, 
        target_path_full, 
        output_path_full,
        output_chimera_path_full,
        thresh_corr=0.9,
        not_found_corr_value=0.0,
        is_log=False,
        log_path=os.getcwd() + os.path.sep + "chimera_logs",
        script_path=os.getcwd() + os.path.sep + "chimera_scripts",
    ):
    """
    Aligns given density map to the target density map: 
    1. Shifts given density map (to_align_path_full) such that its center of mass matches the one 
    of target density map (target_path_full).
    2. Fits shifted density map into the target density map to fully align the maps. 
    3. Checks if the fitting correlation is gerated than the set threshold (thresh_corr). 
    If not, raises RuntimeError.

    Args:
        to_align_path_full - full path to the density we align (including file's name)  
        target_path_full - full path to the target density (including file's name)   
        output_path_full - full path to the output file (including its name) where fitted density will be written
        output_chimera_path_full - full path to the file where Chimera output will be written. 
        NOTE: here it's necessarry to provide this file since fitting scores will be written there
        thresh_corr - threshold value for the correlation between fitted molecule/density and
        target density. If the computed correlation is less than the threshold - raises RunTimeError
        not_found_corr_value - the value which will be returned if the correlation value
        is not found in the Chimera output file
        log_path - path to the folder where the log file will be stored (excluding the file's name which will be created automatically)
        script_path - path to the folder with the python script for Chimera (excluding its name)

    """

    # match center of masses 
    align_center_mass_of_densities(to_align_path_full, target_path_full)

    # fit shifted density map
    with open(output_chimera_path_full, "w") as output_chimera_file: # open file to store fitting scores
        p = fit_into_map_in_chimera(
            to_align_path_full,
            target_path_full,
            output_path_full,
            output_chimera_file,
            is_log=is_log,
            log_path=log_path,
            script_path=script_path
        )
        _, stderr = p.communicate() # catch errors from Chimera subprocess

    if (p.returncode == 0) and (not stderr):
        corr = extract_correlation_from_chimera_file(output_chimera_path_full, not_found_value=not_found_corr_value)
        if corr < thresh_corr:
            raise RuntimeError(f"Error in fitting: computed correlation: {corr} is less than the threshold correlation: {thresh_corr}.")
    else:
        raise RuntimeError(f"Error in fitting: failed to fit through Chimera script. Return code: {p.returncode}, stderr: {stderr}.")


    

def read_density_data_mrc(density_path_full):
    """
    Reads density from the given mrc file, and returns header of the mrc file, density map data
    and voxel size in each direction.

    Args:
        density_path_full - full path to the denstiy mrc file (including its name) 

    Returns:
        density - array with the density
        header - header of the given mrc file
        voxel_size - numpy array containing voxel sizes in all directions (x, y, z)
    """
    assert density_path_full.endswith(".mrc"), "mrc filename must end with .mrc!"

    with mrcfile.open(density_path_full) as mrc:
        density = mrc.data
        header = mrc.header
        voxel_size = np.array(
            [
            mrc.voxel_size.x,
            mrc.voxel_size.y,
            mrc.voxel_size.z,
            ]
        ) 

    return density, header, voxel_size


def calculate_center_of_mass_of_density(density_map, voxel_size=np.ones(3), origin=np.zeros(3)):
    """
    Calculate the center of mass for a 3D cryo-EM density map.
    
    Args:
        density_map - 3D numpy array of the density map
        voxel_size - the size of each voxel (angstroms, for example)
        origin - numpy array with the origin coordinates (in the units of voxel_size)
        
    Returns:
        com - the center of mass (x, y, z) in voxel space or physical space 
        if voxel_size and origin are provided
    """
    # Get the shape of the density map
    shape = density_map.shape
    
    # Create a grid of coordinates for each axis
    x = np.arange(0, shape[0])  # X axis
    y = np.arange(0, shape[1])  # Y axis
    z = np.arange(0, shape[2])  # Z axis
    
    # Calculate the weighted sum for each axis
    total_mass = np.sum(density_map)  # Total mass (sum of all density values)
    
    if total_mass == 0:
        raise ValueError("Density map has no mass (sum of density values is zero).")
    
    x_com = np.sum(density_map * x[:, None, None]) / total_mass
    y_com = np.sum(density_map * y[None, :, None]) / total_mass
    z_com = np.sum(density_map * z[None, None, :]) / total_mass
    
    # The center of mass in voxel coordinates
    com_voxel = np.array([x_com, y_com, z_com])
    
    # Convert to physical coordinates if voxel_size and origin are provided
    com = com_voxel * voxel_size + origin
    
    return com

def calculate_ligand_size(mol):
    """
    Calculate the maximum size (molecular length) of a ligand.
    
    Args:
        mol - RDKit molecule object with 3D coordinates
    
    Returns:
        max_distance - maximum distance between any two atoms (molecular size)
    """
    
    # If the given molecule has more or less than one conformer, it's ambiguous
    # to compute the ligand size.
    assert (
        mol.GetNumConformers() == 1
    ), f"Ambiguity in the ligand size calculation! The RDKit molecule should have only one conformer but has: {mol.GetNumConformers()}"
    
    # Get the 3D coordinates of the atoms
    conformer = mol.GetConformer()
    coords = conformer.GetPositions()
    
    # Calculate the pairwise distances between atoms
    distances = np.linalg.norm(coords[:, None, :] - coords[None, :, :], axis=-1)
    
    # Find the maximum distance (size of the molecule)
    max_distance = np.max(distances)
    
    return max_distance

def log(
    message,
    status="ERROR",
    log_path=os.getcwd() + os.path.sep + "logs",
    log_filename="log.txt",
):
    """
    Creates logs.

    Args:
        message - message to put in the log file
        status - status of the message. Two options: INFO and ERROR
        log_path - path to the directory with log files
        log_filename - name of the log file
    """

    # full path to the log file (including its name)
    path_full = log_path + os.path.sep + log_filename

    # convert message to string in the case if exception was provided
    message = str(message)

    if not message.endswith("\n"):
        message += "\n"

    # log_message includes utc time of the message, status of the message
    # and the message itself
    log_message = (
        datetime.now(timezone.utc).strftime("%H:%M:%S.%f") + " " + status + ": " + message
    )

    with open(path_full, "a") as f:
        f.write(log_message)
        


def create_folder(folder_path):
    """
    Creates a folder if it doesn't exist.

    Args:
        folder_path - path to the folder that shold be created
    """

    if not os.path.exists(folder_path):
        os.mkdir(folder_path)


def split_sdf_file_to_pdbs(
    sdf_path_full, 
    base_pdb_filename, 
    pdb_path=os.getcwd() + os.path.sep + "split_pdb",
    remove_Hs=False,
):
    """
    Splits input .sdf file with some number of molecules to the separate .pdb files
    for each molecule.

    Args:
        sdf_path_full - full path to the input .sdf file
        base_pdb_filename - base filename for output .pdb files
        pdb_path - path to the folder where ouput .pdb files will be stored
        remove_Hs - whether to remove hydrogen atoms from the molecules when reading
    
    Returns:
        n_mols - number of molecules read from the input .sdf file
        pdb_path_full_list - list with full paths to the output pdb files (including their names) 
    """

    assert sdf_path_full.endswith(".sdf"), "sdf filename must end with .sdf!"
    assert base_pdb_filename.endswith(".pdb"), "pdb filename must end with .pdb!"

    mols = []  # list to store RDKit molecule objects

    with Chem.SDMolSupplier(sdf_path_full, sanitize=False, removeHs=remove_Hs) as suppl:
        for mol in suppl:
            if mol is None:
                print(
                    f"Failed to read one of the molecules from the sdf file in: {sdf_path_full}"
                )
                continue
            mols.append(mol)

    pdb_path_full_list = []  # list to store full paths of the output .pdb files   
    for i in range(len(mols)):
        pdb_filename = f"conf_{i + 1}_" + base_pdb_filename
        mol = mols[i]
        pdb_path_full = save_one_conformer_to_pdb(
            mol, -1, pdb_filename, pdb_path=pdb_path
        )
        pdb_path_full_list.append(pdb_path_full)

    n_mols = len(pdb_path_full_list)
    
    return n_mols, pdb_path_full_list


def read_complexes_names_from_file(path_full, skip_header=1, delimiter=","):
    """
    Reads complexes' names from the given file and stores them in numpy array.

    Args:
        path_full - full path to the input file
        skip_header - how many header rows to skip when read the file
        delimiter - delimiter used in the input file (e.g. ',' for .csv)

    Returns:
        complex_names - numpy array with complexes' names
    """

    complex_names = np.genfromtxt(
        path_full, skip_header=skip_header, delimiter=delimiter, dtype=np.str_
    )

    assert (
        len(complex_names.shape) == 1
    ), f"Array with complexes' names should be 1D, but has the shape: {complex_names.shape}"

    return complex_names

def find_file_in_db(complex_name, db_path, filename):
    """
    Tries to find given file for the given complex in the database.
    If the file wasn't found, raises RuntimeError.

    Args:
        complex_name - name (id) of the complex
        db_path - path to the database with complexes files
        filename - name of the file for the given complex
    
    Returns:
        file_path_full - full path to the given file (including its name)
    """

    # path to the folder with the given complex's data
    complex_folder = db_path + os.path.sep + complex_name

    # full path to the file
    file_path_full =  complex_folder + os.path.sep + filename

    # try to find the file 
    if not os.path.isfile(file_path_full):
        raise RuntimeError(f"File: {filename} wasn't found for the complex {complex_name}! Db path: {db_path}.")

    return file_path_full
    

def find_pdb_ligand_file_in_db(complex_name, db_path, base_ligand_name="_ligand.pdb"):
    """
    Finds .pdb ligand file for the given complex in the database.
    If the file wasn't found, raises RuntimeError.

    Args:
        complex_name - name (id) of the complex
        db_path - path to the database with complexes files
        base_ligand_name - base name for the ligand file (used to construct the full name)

    Returns:
        ligand_path_full - full path to the ligand file (icnlduing its name)
    """
    assert base_ligand_name.endswith(".pdb"), "pdb file must end with .pdb!"

    # try to find .pdb file
    ligand_filename = complex_name + base_ligand_name
    ligand_path_full = find_file_in_db(complex_name, db_path, ligand_filename)

    return ligand_path_full

def find_pdb_protein_file_in_db(
    complex_name, db_path, base_protein_name="_protein_processed.pdb"
):
    """
    Finds .pdb protein file for the given complex in the database.
    If the file wasn't found, raises RuntimeError.

    Args:
        complex_name - name (id) of the complex
        db_path - path to the database with complex files
        base_protein_name - base name for the protein file (used to construct the full name)

    Returns:
        full path to the protein file (icnlduing its name)
    """
    assert base_protein_name.endswith(".pdb"), "pdb file must end with .pdb!"

    # try to find .pdb file
    protein_filename = complex_name + base_protein_name
    protein_path_full = find_file_in_db(complex_name, db_path, protein_filename)

    return protein_path_full


def find_mrc_density_file_in_db(
    complex_name, db_path, base_density_name="_ligand.mrc"
):
    """
    Finds .mrc density file for the given complex in the database.
    If the file wasn't found, raises RuntimeError.

    Args:
        complex_name - name (id) of the complex
        db_path - path to the database with complex files
        base_density_name - base name for the density file (used to construct the full name)

    Returns:
        density_path_full - full path to the density file (icnlduing its name)
    """
    assert base_density_name.endswith(".mrc"), "mrc file must end with .mrc!"

    # try to find .mrc file
    density_filename = complex_name + base_density_name
    density_path_full = find_file_in_db(complex_name, db_path, density_filename)

    return density_path_full


def find_txt_file_in_db(
    complex_name, db_path, base_txt_name="_gnina_docked_cc.txt"
):
    """
    Finds a .txt file for the given complex in the database.
    If the file wasn't found, raises RuntimeError.

    Args:
        complex_name - name (id) of the complex
        db_path - path to the database with complex files
        base_txt_name - base name for the txt file (used to construct the full name)

    Returns:
        txt_path_full - full path to the txt file (icnlduing its name)
    """
    assert base_txt_name.endswith(".txt"), "txt file must end with .txt!"

    # try to find .txt file
    txt_filename = complex_name + base_txt_name
    txt_path_full = find_file_in_db(complex_name, db_path, txt_filename)

    return txt_path_full


def apply_args_and_kwargs(fn, args, kwargs):
    """
    Auxiliary function for running functions inside Pool.starmap from multiprocessing.

    Args:
        fn - function to call inside Pool.starmap
        args - positional arguments for fn
        kwargs - keyword arguments for fn

    Returns:
        result of the fn call with args and kwargs
    """
    print(f"Started {fn.__name__} with args: {args} and kwargs: {kwargs}")
    return fn(*args, **kwargs)



