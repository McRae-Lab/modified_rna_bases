🧬 RNA Base Modification Scripts

Herein lay vibe-coded scripts for semi-automating processes in building RNA PDB models — use at your own peril.

These three scripts modify canonical RNA bases in .pdb structures to their modified analogs, preserving geometry.

Each script uses Biopython to modify atom coordinates and naming directly in the PDB file.
They do not rebuild bonds or hydrogens—those are added later by Phenix ReadySet and geometry minimization.

⚙️ Requirements

Python ≥ 3.8

Biopython ≥ 1.8

Phenix ≥ 1.20 (for ReadySet + Minimization)

Install Biopython if needed:

pip install biopython

🚀 Usage

All scripts share the same interface:

python c_to_5mC.py -i input.pdb -o output.pdb
python u_to_pU.py  -i input.pdb -o output.pdb
python u_to_1mpU.py -i input.pdb -o output.pdb


Optional: restrict to specific residues

python c_to_5mC.py -i input.pdb -o output.pdb -r 12,18-21


If -r is omitted, all residues of the target base type (C or U) are modified.

🧩 Step 2 – Add ligand restraints with Phenix ReadySet

After modification, run:

phenix.ready_set modified.pdb


This generates:

output.ligands.cif — custom restraints for modified bases

output.updated.pdb — PDB with added hydrogens and links

🔧 Step 3 – Minimize geometry in Phenix

This step relieves any bond-length or angle outliers after coordinate editing.

phenix.geometry_minimization modified.updated.pdb modified.ligands.cif \
    selection="resname 5MC" max_iterations=300


Examples:

# For 5-methyl-cytidine
phenix.geometry_minimization m5C.updated.pdb m5C.ligands.cif \
    selection="resname 5MC"

# For pseudouridine
phenix.geometry_minimization pU.updated.pdb pU.ligands.cif \
    selection="resname PSU"

# For 1-methyl-pseudouridine
phenix.geometry_minimization 1mpU.updated.pdb 1mpU.ligands.cif \
    selection="resname 1MP"


The resulting minimized model (minimized.pdb) will have:

Correct covalent geometry

Relaxed bond lengths/angles

All hydrogens consistent with the modified bases

🧠 Notes

Residue names follow Phenix/Coot conventions:

5-methyl-cytidine → 5MC

pseudouridine → PSU

1-methyl-pseudouridine → 1MP

Atom names match the Phenix ligand dictionary:

C5 methyl → CM5

N1 methyl → C1

The scripts do not change chain IDs, residue numbers, or occupancy.

Works for RNA fragments, single nucleotides, or entire models.

📜 Citation / Acknowledgment

If you adapt or extend these scripts, acknowledgement would be rad:

McRae E.K.S. (2025).
Center for RNA Therapeutics, Houston Methodist Research Institute.
