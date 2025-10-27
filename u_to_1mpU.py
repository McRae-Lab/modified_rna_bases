#!/usr/bin/env python3
"""
u_to_1mpU.py
------------
Convert uridines (U) to 1-methyl-pseudouridines (B8H).

Steps:
  1) Rigid-body U→Ψ base reorientation (C1′–C5 linkage), same method as u_to_pU.py.
  2) Add N1-methyl carbon (C1) 1.47 Å from N1, in-plane and pointing OUTWARD
     from the ring (robust direction check using the ring centroid).
Conventions:
  - Residue name -> B8H
  - Methyl carbon atom name -> C1

Usage:
  python u_to_1mpU.py -i input.pdb -o output.pdb            # all U residues
  python u_to_1mpU.py -i input.pdb -o output.pdb -r 12

Author: Ewan McRae Lab, 2025
"""

import argparse, math
from Bio.PDB import PDBParser, PDBIO, Atom

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
    if vdot(z,z) < 1e-8:
        tmp = [1,0,0] if abs(x[0]) < 0.9 else [0,1,0]
        z = vnorm(vcross(x, tmp))
    y = vcross(z, x)
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

# ---------- U→Ψ ----------
def convert_U_to_PSU_inplace(residue):
    need = ["N1","C2","C6","C4","C5","C1'"]
    try:
        atoms = {n: residue[n].coord for n in need}
    except KeyError as e:
        print(f"  Skipped {residue.id[1]} (missing {e})")
        return False

    n1, c2, c6, c4, c5, c1p = atoms["N1"], atoms["C2"], atoms["C6"], atoms["C4"], atoms["C5"], atoms["C1'"]
    RU, oU = make_frame(n1, c2, c6)
    RP, oP = make_frame(c5, c4, c6)
    R, t = compose_transform(RU, oU, RP, oP)

    for atom in residue:
        if atom.get_name() in SUGAR:
            continue
        atom.coord = apply_RT(R, t, atom.coord)

    return True  # we'll rename later when methyl added

# ---------- add N1 methyl outward ----------
def add_methyl_N1_outward(residue):
    try:
        n1 = residue["N1"].coord
        c2 = residue["C2"].coord
        c6 = residue["C6"].coord
    except KeyError:
        print(f"  Skipped {residue.id[1]} (missing N1/C2/C6)")
        return False

    # ring centroid to decide outward/inward
    ring_names = ["C2","N3","C4","C5","C6"]
    centroid = [0.0,0.0,0.0]; count=0
    for nm in ring_names:
        if nm in residue:
            for i,x in enumerate(residue[nm].coord):
                centroid[i] += x
            count += 1
    centroid = [x/count for x in centroid] if count>0 else n1

    u = vnorm(vsub(n1, c2))
    v = vnorm(vsub(n1, c6))
    n = vnorm(vcross(u, v))  # plane normal

    def proj(w):
        dot = vdot(w, n)
        wpi = [w[i] - dot * n[i] for i in range(3)]
        L = math.sqrt(vdot(wpi, wpi)) or 1.0
        return [x/L for x in wpi]

    u, v = proj(u), proj(v)
    bis = proj([u[i] + v[i] for i in range(3)])

    # If bisector points toward ring center, flip to ensure "outward"
    to_centroid = vsub(centroid, n1)
    if vdot(bis, to_centroid) > 0:
        bis = [-x for x in bis]

    c1 = vadd(n1, [1.47 * x for x in bis])  # N(sp2)-C(sp3) ~1.47 Å

    residue.add(Atom.Atom(
        "CN1", coord=c1, bfactor=20.0, occupancy=1.0,
        altloc=' ', fullname=" CN1 ", serial_number=0, element="C"
    ))
    return True

# ---------- CLI ----------
def main():
    ap = argparse.ArgumentParser(description="Convert U → B8H (Ψ + N1 methyl).")
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
                    if convert_U_to_PSU_inplace(res) and add_methyl_N1_outward(res):
                        res.resname = "B8H"
                        print(f"  B8H: converted U→Ψ and added C1 at residue {res.id[1]}")
                        n_mod += 1
                elif args.residues and res.id[1] in target_nums:
                    print(f"  Skipped {res.id[1]} (resname {res.get_resname().strip()}, not U)")

    io.set_structure(structure)
    io.save(args.output)
    print(f"\nDone. Modified {n_mod} residues. Wrote: {args.output}")

if __name__ == "__main__":
    main()
