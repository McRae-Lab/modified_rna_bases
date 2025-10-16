#!/usr/bin/env python3
"""
u_to_pU.py
----------
Convert uridines (U) to pseudouridines (PSU) by rigid-body reorienting the base
so the glycosidic linkage changes from C1′–N1 (U) to C1′–C5 (Ψ).

Method:
  1) Build local frames:
       U-frame:  origin=N1, x=N1→C2, plane from (N1,C2,C6)
       Ψ-frame:  origin=C5, x=C5→C4, plane from (C5,C4,C6)
  2) Compute transform mapping Ψ-frame → U-frame; apply to base atoms only.
  3) Rename residue to PSU.

Usage:
  python u_to_pU.py -i input.pdb -o output.pdb            # all U residues
  python u_to_pU.py -i input.pdb -o output.pdb -r 12,33-41

Author: Ewan McRae Lab, 2025
"""

import argparse, math
from Bio.PDB import PDBParser, PDBIO

# ---------- vector helpers ----------
def vsub(a,b): return [a[i]-b[i] for i in range(3)]
def vadd(a,b): return [a[i]+b[i] for i in range(3)]
def vdot(a,b): return a[0]*b[0]+a[1]*b[1]+a[2]*b[2]
def vscale(a,s): return [a[i]*s for i in range(3)]
def vcross(a,b): return [a[1]*b[2]-a[2]*b[1], a[2]*b[0]-a[0]*b[2], a[0]*b[1]-a[1]*b[0]]
def vnorm(a):
    n = math.sqrt(vdot(a,a))
    return a if n == 0 else [a[i]/n for i in range(3)]

def parse_range(r):
    sel=set()
    for t in r.split(','):
        if '-' in t:
            a,b=map(int,t.split('-')); sel.update(range(a,b+1))
        else:
            sel.add(int(t))
    return sel

def make_frame(origin, a, b):
    x = vnorm(vsub(a, origin))
    z = vnorm(vcross(vsub(a, origin), vsub(b, origin)))
    if vdot(z,z) < 1e-8:  # degeneracy guard
        tmp = [1,0,0] if abs(x[0]) < 0.9 else [0,1,0]
        z = vnorm(vcross(x, tmp))
    y = vcross(z, x)
    # columns = x,y,z
    return [[x[0], y[0], z[0]],
            [x[1], y[1], z[1]],
            [x[2], y[2], z[2]]], origin

def matmul3(R, v):
    return [R[0][0]*v[0] + R[0][1]*v[1] + R[0][2]*v[2],
            R[1][0]*v[0] + R[1][1]*v[1] + R[1][2]*v[2],
            R[2][0]*v[0] + R[2][1]*v[1] + R[2][2]*v[2]]

def matT(R):
    return [[R[0][0], R[1][0], R[2][0]],
            [R[0][1], R[1][1], R[2][1]],
            [R[0][2], R[1][2], R[2][2]]]

def matmul3x3(A,B):
    return [[sum(A[i][k]*B[k][j] for k in range(3)) for j in range(3)] for i in range(3)]

def apply_RT(R, t, p): return vadd(matmul3(R, p), t)

def compose_transform(R_to, o_to, R_from, o_from):
    R = matmul3x3(R_to, matT(R_from))
    t = vadd(o_to, matmul3(R, vscale(o_from, -1.0)))
    return R, t

SUGAR = {"P","OP1","OP2","O5'","C5'","C4'","O4'","C3'","O3'","C2'","O2'","C1'"}

# ---------- core op ----------
def convert_U_to_PSU(residue):
    need = ["N1","C2","C6","C4","C5","C1'"]
    try:
        atoms = {n: residue[n].coord for n in need}
    except KeyError as e:
        print(f"  Skipped {residue.id[1]} (missing {e})")
        return False

    n1, c2, c6, c4, c5, c1p = atoms["N1"], atoms["C2"], atoms["C6"], atoms["C4"], atoms["C5"], atoms["C1'"]

    RU, oU = make_frame(n1, c2, c6)  # target = U-frame
    RP, oP = make_frame(c5, c4, c6)  # source = psi-frame
    R, t = compose_transform(RU, oU, RP, oP)

    for atom in residue:
        if atom.get_name() in SUGAR:  # keep sugar fixed
            continue
        atom.coord = apply_RT(R, t, atom.coord)

    residue.resname = "PSU"
    print(f"  pU: converted U→PSU at residue {residue.id[1]}")
    return True

# ---------- CLI ----------
def main():
    ap = argparse.ArgumentParser(description="Convert U → PSU (rigid-body base reorientation).")
    ap.add_argument("-i","--input", required=True)
    ap.add_argument("-o","--output", required=True)
    ap.add_argument("-r","--residues", help="Residue numbers/ranges (e.g., 12,33-41). Default: all U residues.")
    args = ap.parse_args()

    parser = PDBParser(QUIET=True)
    structure = parser.get_structure("rna", args.input)
    io = PDBIO()

    # Selection: all U residues by default
    if args.residues:
        sel = parse_range(args.residues)
        target_nums = sel
        print(f"Targeting residues: {sorted(target_nums)}")
    else:
        target_nums = {res.id[1] for m in structure for ch in m for res in ch if res.get_resname().strip()=="U"}
        print("Targeting ALL uridines (U)")

    n_mod = 0
    for model in structure:
        for chain in model:
            for res in chain:
                if res.get_resname().strip() == "U" and res.id[1] in (target_nums if args.residues else target_nums):
                    if convert_U_to_PSU(res): n_mod += 1
                elif args.residues and res.id[1] in target_nums:
                    print(f"  Skipped {res.id[1]} (resname {res.get_resname().strip()}, not U)")

    io.set_structure(structure)
    io.save(args.output)
    print(f"\nDone. Modified {n_mod} residues. Wrote: {args.output}")

if __name__ == "__main__":
    main()
