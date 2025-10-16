#!/usr/bin/env python3
"""
c_to_5mC.py
-----------
Convert cytidines (C) to 5-methyl-cytidines (5MC) by adding a methyl carbon at C5.

Geometry:
  - Determine the base plane using C4–C5–C6.
  - Place the methyl carbon along the opposite in-plane angle bisector at C5,
    1.50 Å from C5, so angles C4–C5–CM5 and C6–C5–CM5 are ~120°.
Conventions:
  - Residue name -> 5MC
  - Methyl carbon atom name -> CM5 (matches Phenix/Coot dictionaries)

Usage:
  python c_to_5mC.py -i input.pdb -o output.pdb            # all C residues
  python c_to_5mC.py -i input.pdb -o output.pdb -r 16,21-25

Author: Ewan McRae Lab, 2025
"""

import argparse, math
from Bio.PDB import PDBParser, PDBIO, Atom

# ---------- vector helpers ----------
def unit(v):
    n = math.sqrt(sum(x*x for x in v))
    return [x/n for x in v]

def vec(a, b):
    return [b[i]-a[i] for i in range(3)]

def add_vec(a, v, s=1.0):
    return [a[i] + s*v[i] for i in range(3)]

def cross(a,b):
    return [a[1]*b[2]-a[2]*b[1], a[2]*b[0]-a[0]*b[2], a[0]*b[1]-a[1]*b[0]]

def parse_range(rangestr):
    sel=set()
    for token in rangestr.split(','):
        if '-' in token:
            a,b = map(int, token.split('-')); sel.update(range(a,b+1))
        else:
            sel.add(int(token))
    return sel

# ---------- core op ----------
def add_methyl_to_c(residue):
    try:
        c4 = residue["C4"].coord
        c5 = residue["C5"].coord
        c6 = residue["C6"].coord
    except KeyError:
        print(f"  Skipped {residue.id[1]} (missing C4/C5/C6)")
        return False

    # In-plane opposite bisector direction at C5
    u = unit(vec(c5, c4))  # C5->C4
    v = unit(vec(c5, c6))  # C5->C6
    n = unit(cross(u, v))  # plane normal

    def proj_in_plane(w):
        dot = sum(w[i]*n[i] for i in range(3))
        wpi = [w[i]-dot*n[i] for i in range(3)]
        L = math.sqrt(sum(x*x for x in wpi)) or 1.0
        return [x/L for x in wpi]

    u, v = proj_in_plane(u), proj_in_plane(v)
    bis = [u[i] + v[i] for i in range(3)]
    d = proj_in_plane([-bis[i] for i in range(3)])  # opposite bisector

    cm5 = add_vec(c5, d, 1.50)  # C(sp2)-C(sp3) ~1.50 Å

    residue.add(Atom.Atom(
        name="CM5", coord=cm5, bfactor=20.0, occupancy=1.0,
        altloc=' ', fullname=" CM5", serial_number=0, element="C"
    ))
    residue.resname = "5MC"
    print(f"  5mC: added CM5 at residue {residue.id[1]}")
    return True

# ---------- CLI ----------
def main():
    ap = argparse.ArgumentParser(description="Convert C → 5MC (adds CM5).")
    ap.add_argument("-i","--input", required=True)
    ap.add_argument("-o","--output", required=True)
    ap.add_argument("-r","--residues", help="Residue numbers/ranges (e.g., 16,18-22). Default: all C residues.")
    args = ap.parse_args()

    parser = PDBParser(QUIET=True)
    structure = parser.get_structure("rna", args.input)
    io = PDBIO()

    # Selection: all C residues by default
    if args.residues:
        sel = parse_range(args.residues)
        target_nums = sel
        print(f"Targeting residues: {sorted(target_nums)}")
    else:
        target_nums = {res.id[1] for m in structure for ch in m for res in ch if res.get_resname().strip()=="C"}
        print("Targeting ALL cytidines (C)")

    n_mod = 0
    for model in structure:
        for chain in model:
            for res in chain:
                if res.get_resname().strip() == "C" and res.id[1] in (target_nums if args.residues else target_nums):
                    if add_methyl_to_c(res): n_mod += 1
                elif args.residues and res.id[1] in target_nums:
                    print(f"  Skipped {res.id[1]} (resname {res.get_resname().strip()}, not C)")

    io.set_structure(structure)
    io.save(args.output)
    print(f"\nDone. Modified {n_mod} residues. Wrote: {args.output}")

if __name__ == "__main__":
    main()
