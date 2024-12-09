

######### Perform docking ############### 
# from vina import Vina
import os
import shutil
import pandas as pd
import mrcfile
import numpy as np
from subprocess import Popen, PIPE, DEVNULL
from tqdm import tqdm
from rdkit import Chem
from rdkit.Chem import AllChem
from pymol import cmd
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

)


# # Conver PDB to PDBQT files
# def convert_pdb_to_pdbqt(data_dir, data_df):
#     pbar = tqdm(total=len(data_df))
#     # ob_Conversion = openbabel.OBConversion()
#     # ob_Conversion.SetInAndOutFormats("pdb", "pdbqt")
#     for i, row in data_df.iterrows():
#         # cid, pKa = row['pdbid'], float(row['-logKd/Ki'])
#         cid = row['pdbid']
#         protein_path_input = os.path.join(data_dir, cid, f'{cid}_protein_processed.pdb')
#         protein_path_output = os.path.join(data_dir, cid, f'{cid}_protein_processed ')
#         ligand_path_input = os.path.join(data_dir, cid, f'{cid}_ligand.pdb')
#         ligand_path_output = os.path.join(data_dir, cid, f'{cid}_ligand.pdbqt')
   
#         # mol = openbabel.OBMol()
#         # if not ob_Conversion.ReadFile(mol, protein_path_input):
#         #     print(f"Error reading {protein_path_input}")
#         #     return
        
#         # # Write to the output PDBQT file
#         # if not ob_Conversion.WriteFile(mol, protein_path_output):
#         #     print(f"Error writing {protein_path_output}")
#         #     return
        
#         # if not ob_Conversion.ReadFile(mol, ligand_path_input):
#         #     print(f"Error reading {ligand_path_input}")
#         #     return
        
#         # # Write to the output PDBQT file
#         # if not ob_Conversion.WriteFile(mol, ligand_path_output):
#         #     print(f"Error writing {ligand_path_output}")
#         #     return

#         bashCommand_ligand = f"prepare_ligand4 -l {ligand_path_input} -o {ligand_path_output}"   # this needs vina
#         os.system(bashCommand_ligand)

#         bashCommand_receptor = f"prepare_receptor4 -r {protein_path_input} -o {protein_path_output}"
#         os.system(bashCommand_receptor)
    
#         print(f"Conversion successful PROTEIN: {protein_path_input} -> {protein_path_output}")
#         print(f"Conversion successful LIGAND: {ligand_path_input} -> {ligand_path_output}")


# def Vina_docking(data_dir,data_df):
#     pbar = tqdm(total=len(data_df))
#     v = Vina(sf_name='vina')
#     for i, row in data_df.iterrows():
#         # cid, pKa = row['pdbid'], float(row['-logKd/Ki'])
#         cid = row['pdbid']
#         receptor = os.path.join(data_dir, cid, f'{cid}_protein_processed')
#         ligand = os.path.join(data_dir, cid, f'{cid}_ligand')
#         output = os.path.join(data_dir, cid, f'{cid}_vina_out')

#         mrc = mrcfile.open(f'{ligand}.mrc', mode='r+')
#         density_map = mrc.data
#         origin = mrc.header.origin
#         origin = np.array(origin.tolist())

#         # Voxel size in angstroms (modify if known from the MRC header)
#         voxel_size = mrc.voxel_size.x  # Or manually input the voxel size if known

#         center_of_mass = calculate_center_of_mass_of_density(density_map, voxel_size) + origin

#         v.set_receptor(f'{receptor}.pdbqt')
#         v.set_ligand_from_file(f'{ligand}.pdbqt')
        
#         v.compute_vina_maps(center=center_of_mass, box_size=[30, 30, 30])

#         # Score the current pose
#         energy = v.score()
#         print('Score before minimization: %.3f (kcal/mol)' % energy[0])

#         # Minimized locally the current pose
#         energy_minimized = v.optimize()
#         print('Score after minimization : %.3f (kcal/mol)' % energy_minimized[0])
#         v.write_pose(f'{ligand}_minimized.pdbqt', overwrite=True)

#         # Dock the ligand
#         v.dock(exhaustiveness=32, n_poses=20)
#         v.write_poses(f'{output}.pdbqt', n_poses=10, overwrite=True)


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

# def save_frames_from_trajectory_rdkit(data_dir,data_df):
#     pbar = tqdm(total=len(data_df))
#     for i, row in data_df.iterrows():
#         # cid, pKa = row['pdbid'], float(row['-logKd/Ki'])
#         cid = row['pdbid']
#         print("PDB:", cid)
#         # input = os.path.join(data_dir, cid, f'{cid}_gnina_docked.sdf.gz')
#         # bashCommand_ligand = f"gunzip {input}"
#         # os.system(bashCommand_ligand)
        
#         sdf_file = os.path.join(data_dir, cid, f'{cid}_gnina_docked.sdf')
#         output = os.path.join(data_dir, cid, f'{cid}_gnina_docked_frame')

#         suppl = Chem.SDMolSupplier(sdf_file)
#         for idx, mol in enumerate(suppl):
#             if mol is None:
#                 continue
#             # Save the current molecule to a temporary SDF file
#             mol_sdf = f"{output}_{idx}.sdf"
#             writer = Chem.SDWriter(mol_sdf)
#             writer.write(mol)
#             writer.close()

# def save_frames_from_trajectory_pymol(data_dir,data_df):


#     sdf_file_trajectory = os.path.join(data_dir, cid, f'{cid}_gnina_docked.sdf')
#     output = os.path.join(data_dir, cid, f'{cid}_gnina_docked_frame')


#     # Load the structure and trajectory files
#     cmd.load(sdf_file_trajectory, "molecule")  # Load the trajectory into the structure

#     # Get the number of frames in the trajectory
#     num_frames = cmd.count_states("molecule")

#     # Save each frame separately
#     for frame in range(1, num_frames + 1):
#         cmd.frame(frame)
#         output_file = f"{output}_{frame}.mol2"
#         cmd.save(output_file, "molecule")
#         print(f"Saved {output_file}")

#         # Clean up
#     cmd.delete("molecule")

# def run_chimera(path_to_pdb, path_to_density):
#     # NOTE: this function works only for ChimeraX
#     with open('chimera_tmp_script.py', 'w') as f:
#         f.write('''from chimera import runCommand as rc\nrc("open {path1}")\nrc("molmap #0 5")\nrc("open {path2}")\nrc("measure correlation #0.1 #1")\n'''.format(path1 = path_to_pdb, path2 = path_to_density))
#     os.system('chimera --nogui --script chimera_tmp_script.py  > chimera_res.tmp')
#     with open('chimera_res.tmp', 'r') as f2:
#         for line in f2:
#             if "correlation" in line:
#                 cc = float(line.strip().split()[2][:-1])
#                 break
#     return cc


def filter_docked_poses_by_correlation(
    docked_path_full_list,
    target_density_path_full,
    threshold_correlation=0.6,
    not_found_corr_value=0.0,
    chimeraX_output_base_filename="output.txt",
    chimeraX_output_path=os.getcwd() + os.path.sep + "chimeraX_output",
    density_resolution=1.0,
    n_box=16,
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
):
    """
    The main function that generates docked poses, compute cross correlations for them and write the correlations to a file.

    Args:

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



# def get_cross_correlation(data_dir,data_df):
#     pbar = tqdm(total=len(data_df))
#     for i, row in data_df.iterrows():
#         # cid, pKa = row['pdbid'], float(row['-logKd/Ki'])
#         cid = row['pdbid']
#         print("PDB:", cid)
#         path_to_density = os.path.join(data_dir, cid, f'{cid}_ligand.mrc')
#         ligand_docked_files = os.path.join(data_dir, cid, f'{cid}_gnina_docked_frame')
#         CC_output_file_name = os.path.join(data_dir, cid, f'{cid}_gnina_docked_cc.txt')
#         f = open(CC_output_file_name, "w")
#         for idx in range(1,10):
#             try:
#                 pdb_file = f"{ligand_docked_files}_{idx}.mol2"
#                 mrc_file = f"{path_to_density}"
#                 print(pdb_file)
#                 print(mrc_file)
#                 cc = run_chimera(pdb_file, mrc_file)
#                 print(cc)
#                 f.write(str(cc) + "\n")
#             except:
#                 pass
#         f.close()


    
if __name__ == '__main__':
    pass
    # distance = 5
    # input_ligand_format = 'mol2'
    # data_root = '../../../../../PDBbind'
    # data_dir = os.path.join(data_root, 'PDBBind_Zenodo_6408497')
    # # data_df = pd.read_csv(os.path.join(data_root, 'toy_examples.csv'))
    # data_df = pd.read_csv(os.path.join(data_dir, 'PDB_IDs_with_rdkit_length_less_than_16A_succ_gnina.csv'))
    # # data_df = pd.read_csv(os.path.join(data_dir, 'temp.csv'))

    # ## generate pocket within 5 Ångström around ligand 
    # #convert_pdb_to_pdbqt(data_dir,data_df)
    # # Vina_docking(data_dir,data_df)
    # # gnina_docking(data_dir,data_df)
    # # save_frames_from_trajectory_pymol(data_dir,data_df)
    # get_cross_correlation(data_dir,data_df)
    


