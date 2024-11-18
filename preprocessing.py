# %%
import os
import pickle
from rdkit import Chem
from rdkit.Chem import AllChem
import pandas as pd
from tqdm import tqdm
import pymol
from rdkit import RDLogger
RDLogger.DisableLog('rdApp.*')
from openbabel import openbabel
import numpy as np
import mrcfile

# %%
def generate_pocket(data_dir, data_df, distance=5):
    # complex_id = os.listdir(data_dir)
    # for cid in complex_id:
    for i, row in data_df.iterrows():
        cid = row['pdbid']
        complex_dir = os.path.join(data_dir, cid)
        lig_native_path = os.path.join(complex_dir, f"{cid}_ligand.mol2")
        #protein_path= os.path.join(complex_dir, f"{cid}_protein.pdb")
        protein_path= os.path.join(complex_dir, f"{cid}_protein_processed.pdb")

        if os.path.exists(os.path.join(complex_dir, f'Pocket_{distance}A.pdb')):
            continue

        pymol.cmd.load(protein_path)
        pymol.cmd.remove('resn HOH')
        pymol.cmd.load(lig_native_path)
        pymol.cmd.remove('hydrogens')
        pymol.cmd.select('Pocket', f'byres {cid}_ligand around {distance}')
        pymol.cmd.save(os.path.join(complex_dir, f'Pocket_{distance}A.pdb'), 'Pocket')
        pymol.cmd.delete('all')

def generate_complex(data_dir, data_df, distance=5, input_ligand_format='mol2'):
    pbar = tqdm(total=len(data_df))
    for i, row in data_df.iterrows():
        # cid, pKa = row['pdbid'], float(row['-logKd/Ki'])
        cid = row['pdbid']
        complex_dir = os.path.join(data_dir, cid)
        pocket_path = os.path.join(data_dir, cid, f'Pocket_{distance}A.pdb')
        if input_ligand_format != 'pdb':
            ligand_input_path = os.path.join(data_dir, cid, f'{cid}_ligand.{input_ligand_format}')
            ligand_path = ligand_input_path.replace(f".{input_ligand_format}", ".pdb")
            # os.system(f'obabel {ligand_input_path} -O {ligand_path} -d')
            obConversion.ReadFile(mol, ligand_input_path)
            obConversion.WriteFile(mol, ligand_path)
        else:
            ligand_path = os.path.join(data_dir, cid, f'{cid}_ligand.pdb')

        save_path = os.path.join(complex_dir, f"{cid}_{distance}A.rdkit")
        ligand = Chem.MolFromPDBFile(ligand_path, removeHs=True)
        if ligand == None:
            print(f"Unable to process ligand of {cid}")
            continue

        pocket = Chem.MolFromPDBFile(pocket_path, removeHs=True)
        if pocket == None:
            print(f"Unable to process protein of {cid}")
            continue

        complex = (ligand, pocket)
        with open(save_path, 'wb') as f:
            pickle.dump(complex, f)

        pbar.update(1)


# Conver PDB to PDBQT files
def convert_pdb_to_pdbqt(data_dir, data_df):
    pbar = tqdm(total=len(data_df))
    # ob_Conversion = openbabel.OBConversion()
    # ob_Conversion.SetInAndOutFormats("pdb", "pdbqt")
    for i, row in data_df.iterrows():
        # cid, pKa = row['pdbid'], float(row['-logKd/Ki'])
        cid = row['pdbid']
        protein_path_input = os.path.join(data_dir, cid, f'{cid}_protein_processed.pdb')
        protein_path_output = os.path.join(data_dir, cid, f'{cid}_protein_processed.pdbqt')
        ligand_path_input = os.path.join(data_dir, cid, f'{cid}_ligand.pdb')
        ligand_path_output = os.path.join(data_dir, cid, f'{cid}_ligand.pdbqt')
   
        # mol = openbabel.OBMol()
        # if not ob_Conversion.ReadFile(mol, protein_path_input):
        #     print(f"Error reading {protein_path_input}")
        #     return
        
        # # Write to the output PDBQT file
        # if not ob_Conversion.WriteFile(mol, protein_path_output):
        #     print(f"Error writing {protein_path_output}")
        #     return
        
        # if not ob_Conversion.ReadFile(mol, ligand_path_input):
        #     print(f"Error reading {ligand_path_input}")
        #     return
        
        # # Write to the output PDBQT file
        # if not ob_Conversion.WriteFile(mol, ligand_path_output):
        #     print(f"Error writing {ligand_path_output}")
        #     return

        bashCommand_ligand = f"prepare_ligand4 -l {ligand_path_input} -o {ligand_path_output}"
        os.system(bashCommand_ligand)

        bashCommand_receptor = f"prepare_receptor4 -r {protein_path_input} -o {protein_path_output}"
        os.system(bashCommand_receptor)
    
        print(f"Conversion successful PROTEIN: {protein_path_input} -> {protein_path_output}")
        print(f"Conversion successful LIGAND: {ligand_path_input} -> {ligand_path_output}")

obConversion = openbabel.OBConversion()
obConversion.SetInAndOutFormats("mol2", "pdb")
mol = openbabel.OBMol()


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


if __name__ == '__main__':
    distance = 5
    input_ligand_format = 'mol2'
    data_root = '../../../../../PDBbind'
    data_dir = os.path.join(data_root, 'PDBBind_Zenodo_6408497')
    # data_df = pd.read_csv(os.path.join(data_root, 'toy_examples.csv'))
    data_df = pd.read_csv(os.path.join(data_dir, 'PDB_IDs_with_rdkit_length_less_than_16A.csv'))

    ## generate pocket within 5 Ångström around ligand 
    #generate_pocket(data_dir, data_df, distance=distance)
    #generate_complex(data_dir, data_df, distance=distance, input_ligand_format=input_ligand_format)
    # convert_pdb_to_pdbqt(data_dir,data_df)
    # gnina_docking(data_dir,data_df)


# %%
