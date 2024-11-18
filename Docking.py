

######### Perform docking ############### 
# from vina import Vina
import os
from tqdm import tqdm
import pandas as pd
import mrcfile
import numpy as np
from rdkit import Chem
from rdkit.Chem import AllChem
from pymol import cmd

def calculate_center_of_mass_of_density(density_map, voxel_size=1.0):
    """
    Calculate the center of mass for a 3D cryo-EM density map.
    
    Args:
        density_map (numpy.ndarray): 3D numpy array of the density map.
        voxel_size (float): The size of each voxel (angstroms, for example).
        
    Returns:
        numpy.ndarray: The center of mass (x, y, z) in voxel space or physical space if voxel_size is provided.
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
    
    # Convert to physical coordinates if voxel_size is provided
    com_physical = com_voxel * voxel_size
    
    return com_physical



def Vina_docking(data_dir,data_df):
    pbar = tqdm(total=len(data_df))
    v = Vina(sf_name='vina')
    for i, row in data_df.iterrows():
        # cid, pKa = row['pdbid'], float(row['-logKd/Ki'])
        cid = row['pdbid']
        receptor = os.path.join(data_dir, cid, f'{cid}_protein_processed')
        ligand = os.path.join(data_dir, cid, f'{cid}_ligand')
        output = os.path.join(data_dir, cid, f'{cid}_vina_out')

        mrc = mrcfile.open(f'{ligand}.mrc', mode='r+')
        density_map = mrc.data
        origin = mrc.header.origin
        origin = np.array(origin.tolist())

        # Voxel size in angstroms (modify if known from the MRC header)
        voxel_size = mrc.voxel_size.x  # Or manually input the voxel size if known

        center_of_mass = calculate_center_of_mass_of_density(density_map, voxel_size) + origin

        v.set_receptor(f'{receptor}.pdbqt')
        v.set_ligand_from_file(f'{ligand}.pdbqt')
        
        v.compute_vina_maps(center=center_of_mass, box_size=[30, 30, 30])

        # Score the current pose
        energy = v.score()
        print('Score before minimization: %.3f (kcal/mol)' % energy[0])

        # Minimized locally the current pose
        energy_minimized = v.optimize()
        print('Score after minimization : %.3f (kcal/mol)' % energy_minimized[0])
        v.write_pose(f'{ligand}_minimized.pdbqt', overwrite=True)

        # Dock the ligand
        v.dock(exhaustiveness=32, n_poses=20)
        v.write_poses(f'{output}.pdbqt', n_poses=10, overwrite=True)


def calculate_ligand_size(mol):
    """
    Calculate the maximum size (molecular length) of a ligand.
    
    Args:
        mol (rdkit.Chem.Mol): RDKit molecule object with 3D coordinates.
    
    Returns:
        float: Maximum distance between any two atoms (molecular size).
    """
    # Get the 3D coordinates of the atoms
    conformer = mol.GetConformer()
    num_atoms = mol.GetNumAtoms()
    
    # Extract the coordinates of each atom
    coords = np.array([conformer.GetAtomPosition(i) for i in range(num_atoms)])
    
    # Calculate the pairwise distances between atoms
    distances = np.linalg.norm(coords[:, None, :] - coords[None, :, :], axis=-1)
    
    # Find the maximum distance (size of the molecule)
    max_distance = np.max(distances)
    
    return max_distance



def gnina_docking(data_dir,data_df):
    pbar = tqdm(total=len(data_df))
    for i, row in data_df.iterrows():
        # cid, pKa = row['pdbid'], float(row['-logKd/Ki'])
        cid = row['pdbid']
        print("PDB:", cid)
        receptor = os.path.join(data_dir, cid, f'{cid}_protein_processed')
        ligand = os.path.join(data_dir, cid, f'{cid}_ligand')
        output = os.path.join(data_dir, cid, f'{cid}_gnina_docked.sdf.gz')

        mrc = mrcfile.open(f'{ligand}.mrc', mode='r+')
        density_map = mrc.data
        origin = mrc.header.origin
        origin = np.array(origin.tolist())

        # Voxel size in angstroms (modify if known from the MRC header)
        voxel_size = mrc.voxel_size.x  # Or manually input the voxel size if known

        center_of_mass = calculate_center_of_mass_of_density(density_map, voxel_size) + origin

        center_of_mass_x = center_of_mass[0]
        center_of_mass_y = center_of_mass[1]
        center_of_mass_z = center_of_mass[2]

        # PDB files works better than mol2, in mol2 file I was sometimes getting errors 
        # mol = Chem.MolFromMol2File(f'{ligand}.mol2', sanitize=False, removeHs=False)

        try: 
            mol = Chem.MolFromPDBFile(f'{ligand}.pdb', sanitize=False, removeHs=False)
            # Generate 3D coordinates if not already present
            AllChem.EmbedMolecule(mol)
            # Calculate the molecular size
            ligand_size = calculate_ligand_size(mol)
            # print(f"Ligand size: {ligand_size:.2f} Å")
            box_size = ligand_size + 5

            bashCommand_ligand = f"./gnina -r {receptor}.pdbqt -l {ligand}.pdbqt --center_x {center_of_mass_x} --center_y {center_of_mass_y} --center_z {center_of_mass_z} --size_x {box_size} --size_y {box_size} --size_z {box_size} -o {output}"
            os.system(bashCommand_ligand)
        except:
            pass


def save_frames_from_trajectory_rdkit(data_dir,data_df):
    pbar = tqdm(total=len(data_df))
    for i, row in data_df.iterrows():
        # cid, pKa = row['pdbid'], float(row['-logKd/Ki'])
        cid = row['pdbid']
        print("PDB:", cid)
        # input = os.path.join(data_dir, cid, f'{cid}_gnina_docked.sdf.gz')
        # bashCommand_ligand = f"gunzip {input}"
        # os.system(bashCommand_ligand)
        
        sdf_file = os.path.join(data_dir, cid, f'{cid}_gnina_docked.sdf')
        output = os.path.join(data_dir, cid, f'{cid}_gnina_docked_frame')

        suppl = Chem.SDMolSupplier(sdf_file)
        for idx, mol in enumerate(suppl):
            if mol is None:
                continue
            # Save the current molecule to a temporary SDF file
            mol_sdf = f"{output}_{idx}.sdf"
            writer = Chem.SDWriter(mol_sdf)
            writer.write(mol)
            writer.close()


def save_frames_from_trajectory_pymol(data_dir,data_df):
    pbar = tqdm(total=len(data_df))
    for i, row in data_df.iterrows():
        # cid, pKa = row['pdbid'], float(row['-logKd/Ki'])
        cid = row['pdbid']
        print("PDB:", cid)

        # input = os.path.join(data_dir, cid, f'{cid}_gnina_docked.sdf.gz')
        # bashCommand_ligand = f"gunzip {input}"
        # os.system(bashCommand_ligand)

        sdf_file_trajectory = os.path.join(data_dir, cid, f'{cid}_gnina_docked.sdf')
        output = os.path.join(data_dir, cid, f'{cid}_gnina_docked_frame')


        # Load the structure and trajectory files
        cmd.load(sdf_file_trajectory, "molecule")  # Load the trajectory into the structure

        # Get the number of frames in the trajectory
        num_frames = cmd.count_states("molecule")

        # Save each frame separately
        for frame in range(1, num_frames + 1):
            cmd.frame(frame)
            output_file = f"{output}_{frame}.mol2"
            cmd.save(output_file, "molecule")
            print(f"Saved {output_file}")

            # Clean up
        cmd.delete("molecule")

def run_chimera(path_to_pdb, path_to_density):
    with open('chimera_tmp_script.py', 'w') as f:
        f.write('''from chimera import runCommand as rc\nrc("open {path1}")\nrc("molmap #0 5")\nrc("open {path2}")\nrc("measure correlation #0.1 #1")\n'''.format(path1 = path_to_pdb, path2 = path_to_density))
    os.system('chimera --nogui --script chimera_tmp_script.py  > chimera_res.tmp')
    with open('chimera_res.tmp', 'r') as f2:
        for line in f2:
            if "correlation" in line:
                cc = float(line.strip().split()[2][:-1])
                break
    return cc

def get_cross_correlation(data_dir,data_df):
    pbar = tqdm(total=len(data_df))
    for i, row in data_df.iterrows():
        # cid, pKa = row['pdbid'], float(row['-logKd/Ki'])
        cid = row['pdbid']
        print("PDB:", cid)
        path_to_density = os.path.join(data_dir, cid, f'{cid}_ligand.mrc')
        ligand_docked_files = os.path.join(data_dir, cid, f'{cid}_gnina_docked_frame')
        CC_output_file_name = os.path.join(data_dir, cid, f'{cid}_gnina_docked_cc.txt')
        f = open(CC_output_file_name, "w")
        for idx in range(1,10):
            try:
                pdb_file = f"{ligand_docked_files}_{idx}.mol2"
                mrc_file = f"{path_to_density}"
                print(pdb_file)
                print(mrc_file)
                cc = run_chimera(pdb_file, mrc_file)
                print(cc)
                f.write(str(cc) + "\n")
            except:
                pass
        f.close()
    
if __name__ == '__main__':
    distance = 5
    input_ligand_format = 'mol2'
    data_root = '../../../../../PDBbind'
    data_dir = os.path.join(data_root, 'PDBBind_Zenodo_6408497')
    # data_df = pd.read_csv(os.path.join(data_root, 'toy_examples.csv'))
    data_df = pd.read_csv(os.path.join(data_dir, 'PDB_IDs_with_rdkit_length_less_than_16A_succ_gnina.csv'))
    # data_df = pd.read_csv(os.path.join(data_dir, 'temp.csv'))

    ## generate pocket within 5 Ångström around ligand 
    #generate_pocket(data_dir, data_df, distance=distance)
    #generate_complex(data_dir, data_df, distance=distance, input_ligand_format=input_ligand_format)
    #convert_pdb_to_pdbqt(data_dir,data_df)
    # Vina_docking(data_dir,data_df)
    # gnina_docking(data_dir,data_df)
    # save_frames_from_trajectory_pymol(data_dir,data_df)
    get_cross_correlation(data_dir,data_df)
    


