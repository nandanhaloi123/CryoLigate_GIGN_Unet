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
from utils import delete_extension_from_filename, extract_format_from_filename

# %%
# def generate_pocket(data_dir, data_df, distance=5):
#     # complex_id = os.listdir(data_dir)
#     # for cid in complex_id:
#     for i, row in data_df.iterrows():
#         cid = row['pdbid']
#         complex_dir = os.path.join(data_dir, cid)
#         lig_native_path = os.path.join(complex_dir, f"{cid}_ligand.mol2")
#         #protein_path= os.path.join(complex_dir, f"{cid}_protein.pdb")
#         protein_path= os.path.join(complex_dir, f"{cid}_protein_processed.pdb")

#         if os.path.exists(os.path.join(complex_dir, f'Pocket_{distance}A.pdb')):
#             continue

#         pymol.cmd.load(protein_path)
#         pymol.cmd.remove('resn HOH')
#         pymol.cmd.load(lig_native_path)
#         pymol.cmd.remove('hydrogens')
#         pymol.cmd.select('Pocket', f'byres {cid}_ligand around {distance}')
#         pymol.cmd.save(os.path.join(complex_dir, f'Pocket_{distance}A.pdb'), 'Pocket')
#         pymol.cmd.delete('all')

# def generate_complex(data_dir, data_df, distance=5, input_ligand_format='mol2'):
#     pbar = tqdm(total=len(data_df))
#     for i, row in data_df.iterrows():
#         # cid, pKa = row['pdbid'], float(row['-logKd/Ki'])
#         cid = row['pdbid']
#         complex_dir = os.path.join(data_dir, cid)
#         pocket_path = os.path.join(data_dir, cid, f'Pocket_{distance}A.pdb')
#         if input_ligand_format != 'pdb':
#             ligand_input_path = os.path.join(data_dir, cid, f'{cid}_ligand.{input_ligand_format}')
#             ligand_path = ligand_input_path.replace(f".{input_ligand_format}", ".pdb")
#             # os.system(f'obabel {ligand_input_path} -O {ligand_path} -d')
#             obConversion = openbabel.OBConversion()
#             obConversion.SetInAndOutFormats(input_ligand_format, "pdb")
#             mol = openbabel.OBMol()
#             obConversion.ReadFile(mol, ligand_input_path)
#             obConversion.WriteFile(mol, ligand_path)
#         else:
#             ligand_path = os.path.join(data_dir, cid, f'{cid}_ligand.pdb')

#         save_path = os.path.join(complex_dir, f"{cid}_{distance}A.rdkit")
#         ligand = Chem.MolFromPDBFile(ligand_path, removeHs=True)
#         if ligand == None:
#             print(f"Unable to process ligand of {cid}")
#             continue

#         pocket = Chem.MolFromPDBFile(pocket_path, removeHs=True)
#         if pocket == None:
#             print(f"Unable to process protein of {cid}")
#             continue

#         complex = (ligand, pocket)
#         with open(save_path, 'wb') as f:
#             pickle.dump(complex, f)

#         pbar.update(1)



def generate_pocket(data_dir, cid, ligand_filename, protein_filename, distance=5):
    """
    Generates a pocket for the given ligand and protein using Pymol and saves it to a .pdb file.

    Args:
        data_dir - path to the folder with molecules data
        cid - id of the ligand 
        ligand_filename - name of the ligand file (the folder with cid's data can contain multiple ligand files (e.g. from docking))
        protein_filename - name of the protein file (the folder with cid's data can contain multiple protein files)
        distance - threshold distance to generate the pocket

    Returns:
        pocket_path_full - full path to the .pdb file with pocket data
    """
    complex_dir = os.path.join(data_dir, cid)
    lig_native_path = os.path.join(complex_dir, ligand_filename)
    protein_path = os.path.join(complex_dir, protein_filename)

    pocket_path_full = os.path.join(complex_dir, f'Pocket_{distance}A.pdb')

    if not os.path.exists(pocket_path_full):
        pymol.cmd.load(protein_path)
        pymol.cmd.remove('resn HOH')
        pymol.cmd.load(lig_native_path)
        pymol.cmd.remove('hydrogens')
        pymol.cmd.select('Pocket', f'byres {delete_extension_from_filename(ligand_filename)} around {distance}')
        pymol.cmd.save(pocket_path_full, 'Pocket')
        pymol.cmd.delete('all')

    return pocket_path_full

def generate_complex(data_dir, cid, ligand_filename, pocket_path_full, distance=5):
    """
    Generates complex from the given ligand and computed pocket and saves it to a .rdkit file.

    Args:
        data_dir - path to the folder with molecules data
        cid - id of the ligand 
        ligand_filename - name of the ligand file (the folder with cid's data can contain multiple ligand files (e.g. from docking))
        distance - threshold distance to generate the pocket
        pocket_path_full - full path to the file with the computed pocket
    Returns:
        complex_path_full - full path to the .rdkit file with complex data
    """

    complex_dir = os.path.join(data_dir, cid)
    input_ligand_format = extract_format_from_filename(ligand_filename)

    if input_ligand_format != 'pdb':
        ligand_input_path = os.path.join(complex_dir, ligand_filename)
        ligand_path = os.path.join(complex_dir, delete_extension_from_filename(ligand_filename) + ".pdb")
        # os.system(f'obabel {ligand_input_path} -O {ligand_path} -d')
        obConversion = openbabel.OBConversion()
        obConversion.SetInAndOutFormats(input_ligand_format, "pdb")
        mol = openbabel.OBMol()
        obConversion.ReadFile(mol, ligand_input_path)
        obConversion.WriteFile(mol, ligand_path)
    else:
        ligand_path = os.path.join(complex_dir, ligand_filename)

    complex_path_full = os.path.join(complex_dir, f"{cid}_{distance}A.rdkit")
    ligand = Chem.MolFromPDBFile(ligand_path, removeHs=True)
    if ligand == None:
        raise RuntimeError(f"Failed to generate complex: Unable to process ligand of {cid}")

    pocket = Chem.MolFromPDBFile(pocket_path_full, removeHs=True)
    if pocket == None:
        raise RuntimeError(f"Failed to generate complex: Unable to process pocket of {cid}")

    complex = (ligand, pocket)
    with open(complex_path_full, 'wb') as f:
        pickle.dump(complex, f)

    return complex_path_full

def generate_pocket_and_complex(data_dir, cid, ligand_filename, protein_filename, distance=5):
    """
    Combines generate_pocket() and generate_complex() functions.
    Used for the GIGN dataset creation.

    Args:
        data_dir - path to the folder with molecules data
        cid - id of the ligand 
        ligand_filename - name of the ligand file (the folder with cid's data can contain multiple ligand files (e.g. from docking))
        protein_filename - name of the protein file (the folder with cid's data can contain multiple protein files)
        distance - threshold distance to generate the pocket
    Returns:
        complex_path_full - full path to the .rdkit file with complex data
    """
    
    pocket_path_full = generate_pocket(data_dir, cid, ligand_filename, protein_filename, distance=distance)
    complex_path_full = generate_complex(data_dir, cid, ligand_filename, pocket_path_full, distance=distance)

    return complex_path_full

if __name__ == '__main__':
    pass
    # distance = 5
    # input_ligand_format = 'mol2'
    # data_root = '../../../../../PDBbind'
    # data_dir = os.path.join(data_root, 'PDBBind_Zenodo_6408497')
    # data_df = pd.read_csv(os.path.join(data_root, 'toy_examples.csv'))
    # data_df = pd.read_csv(os.path.join(data_dir, 'PDB_IDs_with_rdkit_length_less_than_16A.csv'))

    ## generate pocket within 5 Ångström around ligand 
    #generate_pocket(data_dir, data_df, distance=distance)
    #generate_complex(data_dir, data_df, distance=distance, input_ligand_format=input_ligand_format)


# %%
