import sys
import getopt
import os
from datetime import datetime, timezone

# append the repo's folder to sys for convenient imports
sys.path.append(os.path.dirname(os.path.dirname(os.path.realpath(__file__))))

from chimeraX_scripts.chimeraX_utilities import (
    log,
    delete_extension_from_filename,
    extract_filename_from_full_path,
    run_chimeraX_command as run_com,
)


"""
This script is used to compute density map for the given molecule inside ChimeraX. 
The generated map is cubic with the specified (equaled) number of points in each direction.
The script accepts the following input arguments:
-i: full path to the input molecule file
-r: resolution for the molecule map generation
-b: number of points for the cubic map
-g: grid spacing for the cubic map
-o: full path to the output density file
-l: if provided, contains path to the folder where log file will be written
"""

opts, args = getopt.getopt(sys.argv[1:], "i:r:b:g:o:l:")

molecule_path_full = None  # full path to the input molecule file (including its name)
density_resolution = None  # density map resolution for the molecule (in Angrstrom)
n_box = None  # number of points for the cubic map
grid_spacing = None # grid spacing for the cubic map
density_path_full = None  # full path to the output density file
is_log = False  # whether we should write logs for ChimeraX script
log_path = None  # path to the log folder (excluding log file name)
log_fname = None  # name of the log file
for opt, arg in opts:
    if opt == "-i":
        molecule_path_full = arg
    elif opt == "-r":
        density_resolution = arg
    elif opt == "-b":
        n_box = arg
    elif opt == "-g":
        grid_spacing = arg
    elif opt == "-o":
        density_path_full = arg
    elif opt == "-l":
        is_log = True
        log_path = arg

# check if the arguments are not None after arguments parsing
assert (
    molecule_path_full
), "Path to the input molecule file is None after argument's parsing."

assert density_resolution, "Density resolution is None after argument's parsing."
assert n_box, "Number of points for the cubic map is None after argument's parsing."
assert grid_spacing, "Grid spacing for the cubic map is None after argument's parsing."
assert (
    density_path_full
), "Path to the output density file  is None after argument's parsing."
if is_log:
    assert log_path, "Path to the log folder is None after argument's parsing."

# create log file name if required
if is_log:
    # extract name of the molecule file
    molecule_fname = extract_filename_from_full_path(molecule_path_full)
    log_fname = (
        datetime.now(timezone.utc).strftime("%d_%m_%Y_%H.%M.%S.%f")[:-2]
        + "_density_map_"
        + delete_extension_from_filename(molecule_fname)
        + "_log.txt"
    )

# run ChimeraX commands
if is_log:
    log(
        "Started ChimeraX commands.",
        status="INFO",
        log_path=log_path,
        log_filename=log_fname,
    )

run_com(
    session,
    "open " + molecule_path_full,
    is_log=is_log,
    log_path=log_path,
    log_filename=log_fname,
)  # open molecule file

run_com(
    session,
    "molmap #1 {} cube {} gridSpacing {}".format(density_resolution, n_box, grid_spacing),
    is_log=is_log,
    log_path=log_path,
    log_filename=log_fname,
)  # generate density map for the molecule

run_com(
    session,
    "save {} #2".format(density_path_full),
    is_log=is_log,
    log_path=log_path,
    log_filename=log_fname,
)  # save the map

# if all commands run successfully, log info message and stop ChimeraX
if is_log:
    log(
        "Successfully finished ChimeraX commands! Exiting...",
        status="INFO",
        log_path=log_path,
        log_filename=log_fname,
    )

run_com(
    session,
    "exit",
    is_log=is_log,
    log_path=log_path,
    log_filename=log_fname,
)
