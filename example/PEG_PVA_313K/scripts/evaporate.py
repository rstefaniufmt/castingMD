#!/usr/bin/env python3
"""
Remove solvente acima de um plano.
Água e outros solventes evaporam com base na taxa configurada no solvent.conf
e em sua concentração atual no sistema.
"""

import sys, os, re
import MDAnalysis as mda

if len(sys.argv) != 6:
    print("Uso: evaporate.py in.gro in.top out.gro out.top plane_frac")
    sys.exit(1)

in_gro, in_top, out_gro, out_top, plane_frac = sys.argv[1:]
plane_frac = float(plane_frac)

def load_solvent_rates(config_file="solvent.conf"):
    rates = {}
    if os.path.exists(config_file):
        with open(config_file, "r") as f:
            for line in f:
                line = line.strip()
                if not line or line.startswith("#"): continue
                if "=" in line:
                    key, val = line.split("=", 1)
                    key = key.strip()
                    val = val.split("#", 1)[0].strip()
                    try:
                        rates[key] = float(val)
                    except ValueError:
                        pass
    if not rates:
        rates = {"SOL": 1.0}
    return rates

solvent_rates = load_solvent_rates("solvent.conf")
solvent_names = list(solvent_rates.keys())

u = mda.Universe(in_gro)
Lz = u.dimensions[2]                       # dimensão z da caixa (nm)
cut_sol = plane_frac * Lz                 # limite de evaporação

# ---- coletar informações do sistema ----
system_counts = {}
residues_by_sol = {}
residues_above = []

for resname in solvent_names:
    res = list(u.select_atoms(f"resname {resname}").residues)
    system_counts[resname] = len(res)
    residues_by_sol[resname] = res
    
    # Contar quantos desse solvente estão acima do corte
    above = [r for r in res if r.atoms.center_of_mass(pbc=True)[2] > cut_sol]
    residues_above.extend(above)

n_total_above = len(residues_above)

# ---- calcular as cotas de remoção ----
weights = {}
for resname in solvent_names:
    # Lei de Raoult/Henry: Taxa * Concentração Atual
    weights[resname] = solvent_rates[resname] * system_counts[resname]

total_weight = sum(weights.values())

quotas = {}
if total_weight > 0 and n_total_above > 0:
    exact_amounts = {res: n_total_above * (weights[res] / total_weight) for res in solvent_names}
    
    for res in solvent_names:
        quotas[res] = int(exact_amounts[res])
    
    remainder = n_total_above - sum(quotas.values())
    
    # Distribuir o resto para quem tem a maior parte fracionária
    fractional_parts = {res: exact_amounts[res] - quotas[res] for res in solvent_names}
    for res in sorted(solvent_names, key=lambda x: fractional_parts[x], reverse=True)[:remainder]:
        quotas[res] += 1
else:
    for res in solvent_names:
        quotas[res] = 0

# Garantir que não tentaremos remover mais do que existe
for res in solvent_names:
    quotas[res] = min(quotas[res], system_counts[res])

# ---- selecionar os resíduos a remover ----
rm_residues = []
rm_counts = {res: 0 for res in solvent_names}

for resname in solvent_names:
    if quotas[resname] > 0:
        # Ordenar os resíduos desse solvente por coordenada Z (do maior pro menor)
        res_sorted = sorted(residues_by_sol[resname], key=lambda r: r.atoms.center_of_mass(pbc=True)[2], reverse=True)
        to_remove = res_sorted[:quotas[resname]]
        rm_residues.extend(to_remove)
        rm_counts[resname] = len(to_remove)

keep = u.atoms
if rm_residues:
    # A maneira mais segura e rápida sem estourar limites de string no mda
    # é ir construindo o union de AtomGroups. Mas pra milhares de resíduos, pode ser lento.
    # Usaremos resids, sabendo que num sistema gro normal eles são suficientes.
    # O método map(str) é rápido o suficiente pra alguns milhares.
    rm_ids = [r.resid for r in rm_residues]
    
    # Se houver muitos resíduos, pode dar memory error/parser error no MDAnalysis
    # Uma forma mais robusta é fatiar:
    chunk_size = 5000
    for i in range(0, len(rm_ids), chunk_size):
        chunk_ids = rm_ids[i:i+chunk_size]
        keep = keep.difference(u.select_atoms("resid " + " ".join(map(str, chunk_ids))))

keep.write(out_gro)

# ---- atualizar .top ----
def patch_top(top_in, top_out, new_counts):
    with open(top_in) as fin, open(top_out, "w") as fout:
        for line in fin:
            tokens = line.split()
            # Procurar linhas do tipo 'SOL 3000' no final da topologia
            if tokens and tokens[0] in new_counts:
                if len(tokens) == 2 and tokens[1].isdigit():
                    fout.write(f"{tokens[0]:<10}{new_counts[tokens[0]]:>6}\n")
                    continue
            fout.write(line)

new_counts = {}
for resname in solvent_names:
    new_counts[resname] = system_counts[resname] - rm_counts[resname]

patch_top(in_top, out_top, new_counts)

print(f"Evaporação concluída. Total acima do corte: {n_total_above}")
for resname in solvent_names:
    print(f"Removidos: {rm_counts[resname]} {resname} (Restam: {new_counts[resname]})")

