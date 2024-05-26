from Bio import PDB
from rdkit import Chem
from openbabel import pybel
import os
import random
import string
import warnings


def initialize_parser_io(prot_file):
    """
    Initializes BioPython parser and IO given a PDB protein structure file

    Parameters
    ----------
        prot_file: str
            Path to PDB file without a ligand

    Returns
    -------
        parser: BioPython PDBParser object
        io: BioPython PDB IO object
    """

    if prot_file[-3:].lower() == "pdb":
        parser = PDB.PDBParser()
        io = PDB.PDBIO()
    else:
        raise ValueError("A PDB input protein structure is required")

    return parser, io


def open_sdf(lig_file):
    """
    Opens SDF file

    Parameters
    ----------
        lig_file: str
            Path to ligand SDF, mol2, or PDB file

    Returns
    -------
        mol_objects: List of RDKit Mol object or objects if there are multiple
    """

    mol_objects = []

    for mol in Chem.SDMolSupplier(lig_file):
        if mol is not None:
            mol = Chem.AddHs(mol, addCoords=True)
            mol_objects.append(mol)
        else:
            warnings.warn("Molecule object is NoneType, skipping")

    return mol_objects


def open_pdb(lig_file):
    """
    Opens PDB file

    Parameters
    ----------
        lig_file: str
            Path to ligand SDF, mol2, or PDB file

    Returns
    -------
        mol: RDKit Mol object
    """

    mol = Chem.MolFromPDBFile(lig_file, removeHs=False)

    if mol is not None:
        mol = Chem.AddHs(mol, addCoords=True)
    else:
        warnings.warn("Molecule object is NoneType, skipping")

    return mol


def create_mol(lig_file):
    """
    Creates molecule from ligand input file

    Parameters
    ----------
        lig_file: str
            Path to ligand SDF, mol2, or PDB file

    Returns
    -------
        mol_objects: List of RDKit Mol object or objects if there are multiple
    """

    mol_objects = []

    if lig_file[-3:].lower() == "sdf":
        mol_objects = open_sdf(lig_file=lig_file)

    elif lig_file[-3:].lower() == "pdb":
        mol = open_pdb(lig_file=lig_file)
        mol_objects.append(mol)

    elif lig_file[-4:].lower() == "mol2":
        mol = Chem.MolFromMol2File(lig_file, removeHs=False)

        # If mol is None, try using OpenBabel converter instead to write temporary SDF file
        if mol is None:
            inp = pybel.readfile("mol2", lig_file)
            out_sdf = pybel.Outputfile("sdf", lig_file[:-4] + "sdf", overwrite=True)

            for mol in inp:
                out_sdf.write(mol)

            mol_objects = open_sdf(lig_file=lig_file[:-4] + "sdf")
            # Cleanup
            os.remove(lig_file[:-4] + "sdf")

        else:
            mol_objects.append(mol)

    else:
        raise ValueError("Invalid molecule file type, please use sdf, pdb, or mol2")

    return mol_objects


def create_master_structure(mol, parser, prot_file):
    """
    Combines protein file with RDKit mol object

    Parameters
    ----------
        mol: RDKit Mol object
            Molecule to write to final file
        parser: BioPython Parser object
            Used to read PDB file
        prot_file: str
            Path to PDB file without a ligand

    Returns
    -------
        master_structure: BioPython Structure object of combined ligand and protein
    """

    # Create temporary PDB molecule output with random filename
    tmp_filename = (
        "".join(random.choice(string.ascii_lowercase) for i in range(42)) + ".pdb"
    )
    Chem.MolToPDBFile(mol, tmp_filename)
    pdb_mol_parser = PDB.PDBParser()
    molecule_pdb = pdb_mol_parser.get_structure(file=tmp_filename, id="lig")

    # Create master structure by copying protein structure
    prot_structure = parser.get_structure("prot", prot_file)
    master_structure = prot_structure.copy()

    # Get and set the molecule chain to Z to be added to the master structure object
    chains = list(molecule_pdb.get_chains())
    chains[0].id = "Z"
    chains[0].detach_parent()
    master_structure[0].add(chains[0])

    # Cleanup
    os.remove(tmp_filename)

    return master_structure


def combine(prot_file, lig_file, output_filename=None):
    """
    Combines protein PDB files with ligand SDF, mol2 or PDB files
    The coordinates have to be in the same reference frame for the protein and ligand if a complex is required

    Parameters
    ----------
        prot_file: str
            Path to PDB file without a ligand
        lig_file: str
            Path to ligand SDF, mol2, or PDB file
        output_filename: str
            Optional path to an output file. This has the be the same format as the input protein file (.pdb or .cif)

    Returns
    -------
        None, saves PDB file in the same location with the protein and ligand filenames combined
    """

    prot_file_extension = prot_file[-4:]

    if output_filename is not None:
        base_output_filename = output_filename[:-4]
        prot_file_extension = output_filename[-4:]
        if (
            prot_file_extension[-3:].lower() != "pdb"
        ):
            raise ValueError(
                "Invalid file extension for output file, should be pdb"
            )

    else:
        if lig_file[-4:].lower() == "mol2":
            base_name = os.path.basename(lig_file)
            lig_part = base_name[:-5]
        elif (lig_file[-3:].lower() == "sdf") or (lig_file[-3:].lower() == "pdb"):
            base_name = os.path.basename(lig_file)
            lig_part = base_name[:-4]
        else:
            raise ValueError("Invalid molecule file type, please use sdf, pdb, or mol2")

        base_output_filename = prot_file[:-4] + "_" + lig_part

    # Check if protein is PDB file, initialize parser accordingly
    parser, io = initialize_parser_io(prot_file=prot_file)

    # Get molecule objects
    mol_objects = create_mol(lig_file=lig_file)
    mol_objects = [mol for mol in mol_objects if mol is not None]

    counter = 0
    if len(mol_objects) > 1:
        for mol in mol_objects:
            counter += 1
            master_structure = create_master_structure(
                mol=mol, parser=parser, prot_file=prot_file
            )
            io.set_structure(master_structure)
            io.save(base_output_filename + "_" + str(counter) + prot_file_extension)
    elif len(mol_objects) == 1:
        mol = mol_objects[0]
        master_structure = create_master_structure(
            mol=mol, parser=parser, prot_file=prot_file
        )
        io.set_structure(master_structure)
        io.save(base_output_filename + prot_file_extension)
    else:
        raise ValueError("No molecules could be processed, check inputs")
