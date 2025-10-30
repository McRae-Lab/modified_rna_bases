# RNA Base Modification Scripts 🧬 

Herein lay vibe-coded scripts for semi-automating processes in building RNA PDB models — use at your own peril.

These three scripts modify canonical RNA bases in .pdb structures to their modified analogs, preserving geometry.

Each script uses Biopython to modify atom coordinates and naming directly in the PDB file.
They do not rebuild bonds or hydrogens—those are added later by Phenix ReadySet and geometry minimization.

# Requirements ⚙️ 

Python ≥ 3.8

Biopython ≥ 1.8

Phenix ≥ 1.20 (for ReadySet + Minimization)

Install Biopython if needed:

pip install biopython

# Usage 🚀 

All scripts share the same interface:

python c_to_5mC.py -i input.pdb -o output.pdb
python u_to_pU.py  -i input.pdb -o output.pdb
python u_to_1mpU.py -i input.pdb -o output.pdb


Optional: restrict to specific residues

python c_to_5mC.py -i input.pdb -o 5mC-input.pdb -r 12,18-21


If -r is omitted, all residues of the target base type (C or U) are modified.

# Step 2 – Fix geometry & naming conventions 🔧

For m5C:
Load the model into ChimeraX and initialize isolde, rebuilding residues to match template

Save the pdb - rename the hydrogens to match the cif file with fix_5mC_hnames.sh

For m1pU
Remove the hydrogens
`phenix.reduce m1Y.pdb -Trim > reduced.pdb`

 
# Step 3 – Minimize geometry in Phenix 🔧

This step relieves any bond-length or angle outliers after coordinate editing.

```
phenix.geometry_minimization 5mC-input.pdb 5MC.cif \
    selection="resname 5MC" max_iterations=300
```


Examples:

For 5-methyl-cytidine
```
phenix.geometry_minimization 5mC-input.pdb 5MC.cif \
    selection="resname 5MC max_iterations=300 "
```

For pseudouridine
```
phenix.geometry_minimization pU-input.pdb pU.cif \
    selection="resname PSU max_iterations=300"
```

For 1-methyl-pseudouridine
```
phenix.geometry_minimization reduced.pdb B8H.cif \
    selection="resname B8H max_iterations=300"
```


The resulting minimized model (minimized.pdb) will have:

Correct covalent geometry

Relaxed bond lengths/angles

# Notes 🧠 

Residue names follow Phenix/Coot conventions:

5-methyl-cytidine → 5MC

pseudouridine → PSU

1-methyl-pseudouridine → B8H

Atom names match the Phenix ligand dictionary:

C5 methyl → CM5

N1 methyl → CN1

The scripts do not change chain IDs, residue numbers, or occupancy.

Works for RNA fragments, single nucleotides, or entire models.

# Citation / Acknowledgment 📜 

If you adapt or extend these scripts, acknowledgement would be rad:

McRae E.K.S. (2025).
Center for RNA Therapeutics, Houston Methodist Research Institute.
