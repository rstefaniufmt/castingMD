#!/usr/bin/env python3
"""
Remove solvente acima de um plano.
Água evapora se z > plane*Lz.
Acetona evapora se z > 0.8*plane*Lz (mais rápida).
"""

import sys, re
import MDAnalysis as mda

if len(sys.argv) != 6:
    print("Uso: evaporate.py in.gro in.top out.gro out.top plane_frac")
    sys.exit(1)

in_gro, in_top, out_gro, out_top, plane_frac = sys.argv[1:]
plane_frac = float(plane_frac)

u = mda.Universe(in_gro)
Lz = u.dimensions[2]                       # dimensão z da caixa (nm)
cut_sol = plane_frac * Lz                 # limite para água
cut_ace = 0.8 * plane_frac * Lz           # limite menor p/ acetona

# ---- selecionar resíduos a remover ----
def resid_above(resname, zcut):
    residues = u.select_atoms(f"resname {resname}").residues
    return [r for r in residues if r.atoms.center_of_mass(pbc=True)[2] > zcut]

rm_sol = resid_above("SOL", cut_sol)
rm_ace = resid_above("ACTN", cut_ace)

keep = u.atoms
if rm_sol:
    keep = keep.difference(u.select_atoms("resid " + " ".join(map(str,[r.resid for r in rm_sol]))))
if rm_ace:
    keep = keep.difference(u.select_atoms("resid " + " ".join(map(str,[r.resid for r in rm_ace]))))

keep.write(out_gro)

# ---- atualizar .top (linhas simples tipo 'SOL  3000') ----
def patch_top(top_in, top_out, nsol_out, nace_out):
    with open(top_in) as fin, open(top_out, "w") as fout:
        for line in fin:
            if line.lstrip().startswith("SOL"):
                tokens = line.split()
                fout.write(f"{tokens[0]:<10}{nsol_out:>6}\n")
            elif line.lstrip().startswith("ACTN"):
                tokens = line.split()
                fout.write(f"{tokens[0]:<10}{nace_out:>6}\n")
            else:
                fout.write(line)

nsol_old = len(u.select_atoms("resname SOL").residues)
nace_old = len(u.select_atoms("resname ACTN").residues)
patch_top(in_top, out_top,
          nsol_old - len(rm_sol),
          nace_old - len(rm_ace))

print(f"Removidos: {len(rm_sol)} água; {len(rm_ace)} acetona")

