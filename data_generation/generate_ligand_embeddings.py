import os
import getopt
import sys
import torch
from itertools import repeat
from datetime import datetime, timezone
from multiprocessing.pool import Pool
from rdkit import Chem
from transformers import AutoModel, AutoTokenizer

# append repo path to sys for convenient imports
sys.path.append(os.path.dirname(os.path.dirname(os.path.realpath(__file__))))

from utils import (
    read_molecule,
    find_pdb_ligand_file_in_db,
    log,
    extract_filename_from_full_path,
    delete_extension_from_filename,
    read_complexes_names_from_file,
    create_folder,
    apply_args_and_kwargs,
)


def main(    
    complex_name,
    db_path,
    model,
    tokenizer,
    base_ligand_filename="_ligand.pdb",
    clarifying_embedding_name="_embedding",
    is_main_log=True,
    main_log_filename="log.txt",
    main_log_path=os.path.join(os.getcwd(), "embeddings_main_logs"),
):
    """
    The main function to generate ligand embeddings from SMILES.

    Args:
        complex_name - name (id) of the protein-ligand complex
        db_path - path to the database with the complexes' data
        model - LLM to generate embeddings from SMILES
        tokenizer - required to generate proper input for the model
        base_ligand_filename - base name for the ligand file (used to construct the full name)
        clarifying_embedding_name - additional part of the embedding file name to clarify 
        (needed because we can generate different embeddings for the same complex e.g. with different models)
        is_main_log - whether to write logs for the main function
        main_log_filename - name of the log file for the main function
        main_log_path - path to the folder where main log files will be stored
    """
    try:
        # find ligand file in the database
        ligand_path = find_pdb_ligand_file_in_db(complex_name, db_path, base_ligand_name=base_ligand_filename)

        # read ligand sturcture and convert it to SMILES
        ligand_mol = read_molecule(ligand_path, remove_Hs=True, sanitize=False) # convert ligand file to rdkit molecule
        ligand_mol = Chem.RemoveHs(ligand_mol) # remove hydrogens to get canonical form
        ligand_smiles = Chem.MolToSmiles(ligand_mol) # generate SMILES from the molecule

        # generate embedding for the ligand's SMILES
        inputs = tokenizer(ligand_smiles, padding=True, return_tensors="pt")
        with torch.no_grad():
            outputs = model(**inputs)
            ligand_embedding = outputs.pooler_output

        # save generated embedding to a file
        ligand_fname = extract_filename_from_full_path(ligand_path)
        output_embedding_path = os.path.join(
            db_path, complex_name, f"{delete_extension_from_filename(ligand_fname) + clarifying_embedding_name}.pyg"
        )
        torch.save(ligand_embedding, output_embedding_path)

        if is_main_log:
            log(
                f"Successfully created embeddings {delete_extension_from_filename(ligand_fname) + clarifying_embedding_name}.",
                status="INFO",
                log_path=main_log_path,
                log_filename=main_log_filename,
            )
    except Exception as e:
        if is_main_log:
            log(
                f"Failed to create embeddings {delete_extension_from_filename(ligand_fname) + clarifying_embedding_name}: {e}",
                status="ERROR",
                log_path=main_log_path,
                log_filename=main_log_filename,
            )

if __name__ == "__main__":

    db_path = '/mnt/cephfs/projects/CryoLigate/PDBbind/PDBBind_Zenodo_6408497' # path to the database

    # load molecule names
    complex_names_csv = (
        db_path + os.path.sep + "PDB_IDs_with_rdkit_length_less_than_24A_nbonds_less_than_10.csv"
    )
    complex_names = read_complexes_names_from_file(complex_names_csv)
    n_complexes = len(complex_names)
    print(f"Computing embeddings for {n_complexes} complexes....")
    
    # Load model and tokenizer
    model_name = "ibm/MoLFormer-XL-both-10pct"
    tokenizer = AutoTokenizer.from_pretrained(model_name, trust_remote_code=True)
    model = AutoModel.from_pretrained(model_name, trust_remote_code=True)

    # create log folder
    is_main_log = True
    main_log_path = None
    if is_main_log:
        main_log_filename = (
            datetime.now(timezone.utc).strftime("%d_%m_%Y_%H.%M.%S.%f")[:-2]
            + "_main_log.txt"
        )
        main_log_path = os.path.join(os.getcwd(), "embeddings_main_logs")
        create_folder(main_log_path)

    # specify keyword arguments for the main function
    base_ligand_filename = "_ligand.pdb"
    clarifying_embedding_name = "_embedding"

    for complex_name in complex_names:
        main(
           complex_name,
           db_path,
           model,
           tokenizer,
           base_ligand_filename=base_ligand_filename,
           clarifying_embedding_name=clarifying_embedding_name,
           is_main_log=is_main_log,
           main_log_filename=main_log_filename,
           main_log_path=main_log_path
        )