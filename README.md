# prot-lig-combine
Combines given protein and ligand into a PDB file

## Set up
Create a conda environment and install the packages below:

```bash
conda install conda-forge::rdkit
conda install conda-forge::biopython
conda install conda-forge::openbabel
```

Install the module using pip:

```bash
pip install .
```

## Usage
To create a combined file see the example below which creates an output PDB file in the input location:

```python
from prot_lig_combine import combine

combine.combineProtLig(prot_file=<pdb_or_mmcif_file>, lig_file=<mol2_sdf_or_pdb_molecule_file>)
```
