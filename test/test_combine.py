from prot_lig_combine import combine

import os
from rdkit import Chem
from Bio import PDB

def test_open_sdf():
    mols = combine.open_sdf("1stp_lig.sdf")

    assert len(mols) == 1
    assert type(mols[0]) == Chem.rdchem.Mol

def test_open_pdb():
    mol = combine.open_pdb("1stp_lig.pdb")

    assert type(mol) == Chem.rdchem.Mol

def test_create_mol():
    mols_sdf = combine.create_mol("1stp_lig.sdf")
    mols_pdb = combine.create_mol("1stp_lig.pdb")
    mols_mol2 = combine.create_mol("1stp_lig.mol2")

    assert type(mols_sdf[0]) == Chem.rdchem.Mol
    assert type(mols_pdb[0]) == Chem.rdchem.Mol
    assert type(mols_mol2[0]) == Chem.rdchem.Mol

def test_create_master_structure_pdb():
    parser, io = combine.initialize_parser_io(prot_file="1stp_prot.pdb")

    print(type(parser))
    
    assert type(parser) == PDB.PDBParser and type(io) == PDB.PDBIO

def test_combine():
    combine.combine(prot_file="1stp_prot.pdb", lig_file="1stp_lig.sdf", output_filename="1stp_comb_sdf.pdb")
    combine.combine(prot_file="1stp_prot.pdb", lig_file="1stp_lig.pdb", output_filename="1stp_comb_pdb.pdb")
    combine.combine(prot_file="1stp_prot.pdb", lig_file="1stp_lig.mol2", output_filename="1stp_comb_mol2.pdb")

    # Check if file exists
    assert os.path.isfile("1stp_comb_sdf.pdb") and os.path.isfile("1stp_comb_pdb.pdb") and os.path.isfile("1stp_comb_mol2.pdb")

    # Open PDB files and check that they contain chain Z (ligand)
    parser = PDB.PDBParser()
    struct_sdf = parser.get_structure(file="1stp_comb_sdf.pdb", id="1stp_sdf")
    struct_pdb = parser.get_structure(file="1stp_comb_pdb.pdb", id="1stp_pdb")
    struct_mol2 = parser.get_structure(file="1stp_comb_mol2.pdb", id="1stp_mol2")

    assert "Z" in [i.id for i in struct_sdf.get_chains()] and "Z" in [i.id for i in struct_pdb.get_chains()] and "Z" in [i.id for i in struct_mol2.get_chains()]

    # Remove files if exist
    os.remove("1stp_comb_sdf.pdb")
    os.remove("1stp_comb_pdb.pdb")
    os.remove("1stp_comb_mol2.pdb")