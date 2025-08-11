# Casting MD: Molecular Dynamics Simulation of Polymer Film Formation via Solvent Evaporation

This repository contains scripts and example files for performing molecular dynamics (MD) simulations of polymer films using GROMACS. The main focus is to model the **casting process**, i.e., the progressive evaporation of solvent in systems containing polymers such as Carbohydrates, PVA, PEG, and bioactive compounds like quercetin.

## 🔧 Project Structure

```
casting_md/
│
├── base_scripts/             # Main reusable scripts for evaporation
│   ├── evaporate.py          # Python script that removes solvent above a cutting plane
│   └── run_evaporate_loop.sh  # Bash script that performs multiple evaporation cycles
│
├── example/                  # Usage examples with different systems
│   ├── PEG_PVA_313K/         # Simulation with PEG and PVA at 313 K
│   └── PVA_quercetin_308K/   # Simulation with PVA + quercetin at 308 K
│
├── mdp/                      # .mdp parameter files organized by stage
│   ├── initial/              # Minimization, NVT, NPT initial stages
│   └── evaporate/            # Parameters for the evaporation cycles
│
└── tutorial.pdf              # Basic tutorial (under construction or for reference)
```

## 🚀 Running the Evaporation Simulation

The evaporation is done in cycles, removing solvent molecules (e.g., water, ethanol or acetone) that are **above a cutting plane along the Z axis** of the box. This mimics the progressive drying of the system.

### 1. `evaporate.py` — Solvent removal script

This Python script uses **MDAnalysis** to identify and remove solvent residues located above a threshold `z` value (defined as a fraction of the box height).

**Usage:**
```bash
python3 evaporate.py in.gro in.top out.gro out.top PLANE
```

- `in.gro` — current structure
- `in.top` — original topology
- `out.gro` — new structure without solvent
- `out.top` — new topology with updated solvent count
- `PLANE` — fraction of the Z-axis above which solvent will be removed (e.g., 0.9 means z > 90%)

### 2. `run_evaporate_loop.sh` — Automated evaporation loop

This script runs a **cycle-based simulation**, where in each cycle:
- A short MD simulation is executed (e.g., 1-5 ns). Please edit md.mdp file to set the desired MD simulation time (default is 1ns)
- Solvent is removed above the `z` plane
- Structure/topology is updated for the next cycle
- A short re-minimization and re-equilibration (100 ps NVT and 250-500 ps NPT) is run for the next cycle. Please, edit nvt_eq.mdp and npt_eq.mdp configuration files to set the desired time for temperature and pressure re-equilibration for your system (default is 100 ps and 500 ps for NVT and NPT, respectively). 

**Main settings:**
```bash
SEG_NS=1            # duration of each cycle in ns
CYCLES=20           # total number of evaporation cycles
MDP=md.mdp          # simulation parameters
TOP=film.top        # initial topology
GRO=npt.gro         # initial structure
PYTHON=python3      # Python with MDAnalysis installed
plane=0.95          # height fraction of box for removal
```

**Execution:**
```bash
bash run_evaporate_loop.sh
```

At the end, you will have multiple structures `stepN_last.gro` with varying levels of evaporation, suitable for final analysis of the film structure. Moreover, you will have final_film.gro after final anneling. 

## 📁 Directory Usage

- `example/PEG_PVA_313K/`: contains `.gro`, `.itp`, `.mdp` and topology files to run the PEG+PVA system.
- `example/PVA_quercetin_308K/`: variation using PVA and quercetin.
- `mdp/`: includes ready-to-use minimization, NVT/NPT and evaporation MD parameter files.
- `gromos54a7_atb.ff/`: GROMOS 54a7 force field converted from ATB, used for polymers and solvents.
- `scripts/` (in some examples): contains system-specific versions of `evaporate.py`.

## 📌 Requirements

- [GROMACS](https://manual.gromacs.org/)
- [MDAnalysis](https://www.mdanalysis.org/) (Python)
- Bash Shell
- Python 3.6+

## 📄 Reference

This script collection was developed for simulating polymer film formation via solvent evaporation, enabling controlled drying within molecular dynamics simulations.

If you use this repository or its content, please cite:

Stefani, R., & Luz, R. (2025). CastingMD: Solvent evaporation molecular dynamics for polymeric film formation (Version 1.0). https://doi.org/10.5281/zenodo.16728915

@misc{stefani2025castingmd,
  author       = {Stefani, Ricardo and Luz, Rudiere},
  title        = {CastingMD: Solvent evaporation molecular dynamics for polymeric film formation},
  year         = {2025},
  publisher    = {Zenodo},
  version      = {1.0},
  doi          = {10.5281/zenodo.16728915},
  url          = {https://doi.org/10.5281/zenodo.16728915 }
  }


---
