# CastingMD: Molecular Dynamics Simulation of Polymer Film Formation via Solvent Evaporation

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.16728915.svg)](https://doi.org/10.5281/zenodo.16728915)

This repository contains automated scripts for performing molecular dynamics (MD) simulations of polymer film formation using GROMACS. The workflow models the **casting process** through progressive solvent evaporation in systems containing polymers (carbohydrates, PVA, PEG) and bioactive compounds (e.g., quercetin).

## 🔧 Project Structure

```
casting_md/
│
├── base_scripts/              # Main reusable scripts for evaporation
│   ├── evaporate.conf         # Configuration file
│   ├── evaporate.py           # Python script for solvent removal
│   ├── evaporate_optimized.sh # Automated workflow (modular version)
│   ├── prepare_evaporation.py # MDP file generator
│   └── run_evaporate_loop_final.sh # Legacy loop script
│
├── example/                   # Usage examples with different systems
│   ├── chitosan/              # Simulation with Chitosan
│   ├── PEG_PVA_313K/          # Simulation with PEG and PVA at 313 K
│   ├── PVA_quercetin__casting_308K/ # Simulation with PVA + quercetin at 308 K
│   └── PVA_Quercetin_processed/ # Processed PVA + quercetin files
│
├── mdp/                       # .mdp parameter files (legacy templates)
│   ├── initial/               # Minimization, NVT, NPT initial stages
│   └── evaporate/             # Parameters for evaporation cycles
│
├── CHANGELOG.md               # Version history and bug fixes
├── LICENSE                    # Project license
├── README.md                  # Main documentation
└── tutorial.pdf               # Basic tutorial
```

## 📋 Overview

The evaporation workflow simulates film casting through iterative MD cycles where solvent molecules are progressively removed, mimicking the drying process. Each cycle consists of:

1. Production MD simulation
2. Solvent evaporation (removal above Z-axis cutting plane)
3. System re-equilibration (EM → NVT → NPT)
4. Box height adjustment for stability

## 🚀 Quick Start

### Prerequisites

**Required Software:**
- [GROMACS](https://manual.gromacs.org/) ≥ 2020.1
- [MDAnalysis](https://www.mdanalysis.org/) (Python library)
- Python 3.6+
- Bash shell
- Standard Unix tools (awk, bc)

**Required Input Files:**
1. Equilibrated system from prior NPT simulation:
   - `npt.gro` - Structure file
   - `npt.cpt` - Checkpoint file
2. Topology file: `film.top`
3. Configuration file: `evaporate.conf`

### Step-by-Step Execution

#### 1. Configure Simulation Parameters

Edit `evaporate.conf` with your system parameters:

```bash
# Simulation time parameters
NVT_PS=100        # NVT equilibration (100 ps)
NPT_PS=1000       # NPT equilibration (1 ns)
SEG_NS=1          # Production MD per cycle (1 ns)
ANNEAL_NS=10      # Final annealing (10 ns)
CYCLES=20         # Number of evaporation cycles

# Thermodynamic parameters
TEMP=308          # Simulation temperature (K)
PRESSURE=1.0      # Simulation pressure (bar)
ANNEAL_TEMP=298   # Final annealing temperature (K)

# MD integration parameters
DT=0.001          # Time step (ps) - 1 fs
RLIST=1.2         # Neighbor list cutoff (nm)
RVDW=1.2          # Van der Waals cutoff (nm)
RCOULOMB=1.2      # Coulomb cutoff (nm)

# Input files
TOP=film.top      # Your topology file
GRO=npt.gro       # Initial structure
CPT=npt.cpt       # Initial checkpoint
```

#### 2. Generate MDP Files

Run the preparation script to create all necessary MDP parameter files:

```bash
python3 prepare_evaporation.py
```

This automatically generates:
- ✓ `em.mdp` - Energy minimization (50,000 steps max)
- ✓ `nvt_eq.mdp` - NVT equilibration (NVT_PS duration)
- ✓ `npt_eq.mdp` - NPT equilibration (NPT_PS duration)
- ✓ `md.mdp` - Production MD (SEG_NS × 1,000,000 steps)
- ✓ `md_anneling.mdp` - Final annealing (TEMP → ANNEAL_TEMP)

#### 3. Run Simulation

Execute the automated workflow:

```bash
bash evaporate_optimized.sh
```

The script will:
- Check for required MDP files
- Execute all evaporation cycles
- Perform final annealing
- Generate complete log file

## 📊 Configuration Parameters

### Simulation Time Control

| Parameter | Default | Description |
|-----------|---------|-------------|
| `NVT_PS` | 100 | NVT equilibration duration (picoseconds) |
| `NPT_PS` | 1000 | NPT equilibration duration (picoseconds) |
| `SEG_NS` | 1 | Production MD duration per cycle (nanoseconds) |
| `ANNEAL_NS` | 10 | Final annealing duration (nanoseconds) |
| `CYCLES` | 20 | Total number of evaporation cycles |

### Thermodynamic Parameters

| Parameter | Default | Description |
|-----------|---------|-------------|
| `TEMP` | 308 | Simulation temperature (Kelvin) |
| `PRESSURE` | 1.0 | Simulation pressure (bar) |
| `ANNEAL_TEMP` | 298 | Final annealing target temperature (Kelvin) |

### MD Integration Parameters

| Parameter | Default | Description |
|-----------|---------|-------------|
| `DT` | 0.001 | Integration time step (ps = 1 fs) |
| `RLIST` | 1.2 | Neighbor list cutoff (nm) |
| `RVDW` | 1.2 | Van der Waals cutoff (nm) |
| `RCOULOMB` | 1.2 | Coulomb cutoff (nm) |

### Output Frequencies

| Parameter | Default | Description |
|-----------|---------|-------------|
| `NVT_OUT_FREQ` | 5000 | NVT output frequency (steps, ~5 ps for dt=0.001) |
| `NPT_OUT_FREQ` | 500 | NPT output frequency (steps, ~0.5 ps for dt=0.001) |
| `MD_OUT_FREQ` | 100000 | Production MD output frequency (steps, ~100 ps for dt=0.001) |

### File Paths

| Parameter | Default | Description |
|-----------|---------|-------------|
| `TOP` | film.top | Topology file |
| `GRO` | npt.gro | Initial structure |
| `CPT` | npt.cpt | Initial checkpoint |
| `LOG_FILE` | evaporate.log | Simulation log |
| `GMX_LOG_FILE` | gromacs_output.log | Log file for GROMACS commands |
| `OUT_FILM` | final_film | Final output name |
| `PYTHON` | python3 | Python executable |

## 🔄 Workflow Details

### Cycle Structure

Each evaporation cycle performs the following steps:

#### 1. Production MD
- Duration: `SEG_NS` nanoseconds
- Uses semi-isotropic pressure coupling
- Outputs every `MD_OUT_FREQ` steps (default: ~100 ps)
- Extracts final frame for evaporation

#### 2. Solvent Evaporation
- Calls `evaporate.py` to remove solvent molecules
- Removes residues above Z-axis cutting plane
- Cutting plane fraction: `(0.01)^(1/CYCLES)`
- Updates topology file automatically

#### 3. Box Height Check
- Ensures Z-dimension ≥ 2×rlist + 0.3 nm
- Prevents "box too small" simulation errors
- Automatically adjusts if necessary
- Maintains system consistency

#### 4. System Re-equilibration
- **Energy Minimization** (50,000 steps maximum)
- **NVT Equilibration** (`NVT_PS` picoseconds, default 100 ps)
- **NPT Equilibration** (`NPT_PS` picoseconds, default 1 ns)

### Final Annealing

After all evaporation cycles:
- Cools system from `TEMP` to `ANNEAL_TEMP`
- Duration: `ANNEAL_NS` nanoseconds (default 10 ns)
- Uses simulated annealing protocol
- Produces final relaxed film structure

## 📁 Output Files

### Per Cycle (i = 0 to CYCLES-1)

```
step{i}.tpr/gro/edr/log/cpt    # Production MD files
step{i}_last.gro               # Final frame from MD
step{i}_post.gro/top           # After solvent evaporation
step{i}_min.gro                # After energy minimization
step{i}_nvt.gro/cpt            # After NVT equilibration
step{i}_npt.gro/cpt            # After NPT equilibration
box_fixed_{i}.gro              # If box height was corrected
```

### Final Output

```
final_film.tpr/gro/edr/log/cpt # Annealed final film
evaporate.log                   # Complete workflow log
```

## 🛠️ Core Scripts

### 1. `evaporate.py` — Solvent Removal Script

Uses **MDAnalysis** to identify and remove solvent residues above a Z-axis threshold. Supports complex solvent mixtures (binary/ternary) by reading specific evaporation rates from an optional `solvent.conf` file (e.g., `SOL = 1.0`, `ACTN = 1.4`).

**Usage:**
```bash
python3 evaporate.py in.gro in.top out.gro out.top PLANE
```

**Arguments:**
- `in.gro` — Current structure
- `in.top` — Original topology
- `out.gro` — Output structure without solvent
- `out.top` — Updated topology with corrected solvent count
- `PLANE` — Fraction of Z-axis above which solvent is removed (e.g., 0.9 = remove above 90% height)

### 2. `prepare_evaporation.py` — MDP File Generator

Generates all GROMACS parameter files based on `evaporate.conf` settings.

**Features:**
- Reads configuration from `evaporate.conf`
- Calculates `nsteps` automatically:
  - NVT: `NVT_PS / DT`
  - NPT: `NPT_PS / DT`
  - MD: `SEG_NS × 1,000,000`
  - Annealing: `ANNEAL_NS × 1,000,000`
- Generates 5 MDP files with correct parameters
- Creates example config if none exists

**Usage:**
```bash
python3 prepare_evaporation.py
```

### 3. `evaporate_optimized.sh` — Main Workflow Script

Automated, modular bash script that orchestrates the entire simulation.

**Key Features:**
- Modular function-based design
- Automatic cycle management
- Box height safety checks
- Comprehensive error handling
- Resume capability (skips completed cycles)
- Detailed logging

**Usage:**
```bash
bash evaporate_optimized.sh
```

## 📝 Examples

### Example 1: Standard Film Casting (PVA + Quercetin at 308 K)

```bash
# Configure
cat > evaporate.conf << EOF
SEG_NS=1
CYCLES=20
TEMP=308
PRESSURE=1.0
TOP=pva_quercetin.top
GRO=npt_equilibrated.gro
CPT=npt_equilibrated.cpt
EOF

# Generate MDP files
python3 prepare_evaporation.py

# Run simulation
bash evaporate_optimized.sh
```

### Example 2: PEG + PVA at 313 K with Extended Equilibration

```bash
# Configure with longer equilibration
cat > evaporate.conf << EOF
NVT_PS=200
NPT_PS=2000
SEG_NS=2
CYCLES=15
TEMP=313
PRESSURE=1.0
TOP=peg_pva.top
GRO=npt.gro
CPT=npt.cpt
EOF

# Generate and run
python3 prepare_evaporation.py
bash evaporate_optimized.sh
```

### Example 3: High-Temperature Rapid Evaporation

```bash
cat > evaporate.conf << EOF
SEG_NS=5
CYCLES=10
TEMP=350
ANNEAL_TEMP=273
ANNEAL_NS=20
MD_OUT_FREQ=50000  # More frequent output
EOF

python3 prepare_evaporation.py
bash evaporate_optimized.sh
```

## 🔧 Customization

### Adjusting Simulation Times

Simply edit `evaporate.conf` and regenerate MDP files:

```bash
# Modify parameters
vim evaporate.conf

# Regenerate MDP files with new settings
python3 prepare_evaporation.py

# The script automatically updates:
# - md.mdp with new nsteps
# - nvt_eq.mdp with new duration
# - npt_eq.mdp with new duration
# - md_anneling.mdp with new cooling schedule
```

### Changing Evaporation Rate

Edit the evaporation plane calculation in `evaporate_optimized.sh`:

```bash
# More aggressive evaporation (faster drying)
readonly EVAP_PLANE=$(awk -v n="${CYCLES}" 'BEGIN {printf "%.6f", (0.001)^(1/n)}')

# Less aggressive evaporation (slower drying)
readonly EVAP_PLANE=$(awk -v n="${CYCLES}" 'BEGIN {printf "%.6f", (0.1)^(1/n)}')
```

### Manual MDP Editing

You can directly edit generated MDP files for fine-tuning:

```bash
# After generating files
python3 prepare_evaporation.py

# Edit specific parameters
vim md.mdp              # Adjust production MD settings
vim nvt_eq.mdp          # Modify NVT coupling parameters
vim npt_eq.mdp          # Change pressure coupling
vim md_anneling.mdp     # Customize cooling schedule
```

## 🐛 Troubleshooting

### Error: "Missing MDP files"

**Solution:** Run the MDP generator first:
```bash
python3 prepare_evaporation.py
```

### Error: Simulation crashes with "box too small"

**Solution:** The script should handle this automatically, but if crashes persist:
- Increase `RLIST_BUFFER` in `evaporate_optimized.sh` (default: 0.3 nm)
- Check for extremely aggressive evaporation settings
- Verify initial box dimensions are adequate

### Error: Checkpoint not found

**Solution:**
- Ensure previous cycle completed successfully
- Check `evaporate.log` for errors
- Verify file existence: `step{i}_npt.cpt`
- May need to restart from last successful cycle

### Error: Topology mismatch

**Solution:**
- Ensure `evaporate.py` correctly updates topology
- Verify number of atoms in `.gro` matches topology
- Check SOL count is updated properly
- Review consistency check messages in log

### Simulation seems stuck

**Solution:**
- Check GROMACS output for warnings
- Monitor `evaporate.log` for progress
- Verify system hasn't run out of solvent
- Consider adjusting equilibration times

## 🔄 Resuming Simulations

The workflow features a **Smart Resume Capability**. If a simulation is abruptly interrupted (e.g., power outage or manual cancellation), the script will automatically detect existing checkpoint files (`.cpt`) and resume the exact step (Minimization, NVT, NPT, or MD) using `gmx mdrun -cpi -append`. It will skip already completed cycles and substeps seamlessly.

```bash
# If simulation was interrupted
# Simply re-run the script - it will pick up exactly where it left off!
bash evaporate_optimized.sh
```

To restart from scratch:
```bash
# Remove all generated files from previous runs
rm -f step* box_fixed* final_film* gromacs_output.log evaporate.log
bash evaporate_optimized.sh
```

## 📚 Directory Usage

- **`example/PEG_PVA_313K/`**: Complete system setup for PEG+PVA at 313 K
  - Includes `.gro`, `.itp`, `.top` files
  - Ready-to-use example configuration
  
- **`example/PVA_quercetin_308K/`**: PVA with quercetin bioactive compound
  - Shows bioactive film formation
  - Different temperature protocol

- **`mdp/`**: Legacy MDP template files
  - `initial/`: EM, NVT, NPT for initial setup
  - `evaporate/`: Parameters for evaporation cycles
  - **Note:** New workflow uses auto-generated MDP files

- **`gromos54a7_atb.ff/`**: GROMOS 54a7 force field
  - Converted from ATB database
  - Used for polymers and solvents
  - Compatible with bioactive compounds

## 📖 Citation

If you use this workflow in your research, please cite:

**Stefani, R., & Luz, R. (2025).** *CastingMD: Solvent evaporation molecular dynamics for polymeric film formation* (Version 1.0). Zenodo. https://doi.org/10.5281/zenodo.16728915

### BibTeX

```bibtex
@misc{stefani2025castingmd,
  author       = {Stefani, Ricardo and Luz, Rudiere},
  title        = {CastingMD: Solvent evaporation molecular dynamics for 
                  polymeric film formation},
  year         = {2025},
  publisher    = {Zenodo},
  version      = {1.0},
  doi          = {10.5281/zenodo.16728915},
  url          = {https://doi.org/10.5281/zenodo.16728915}
}
```

## 📄 License

This project is open-source. Please refer to the LICENSE file for details.

## 🤝 Contributing

Contributions are welcome! Please:
1. Fork the repository
2. Create a feature branch
3. Submit a pull request with detailed description

## 📧 Contact

For questions, issues, or suggestions:
- Open an issue on GitHub
- Check existing documentation
- Review example directories

## 🎯 Version History

### Version 2.1 (Current)
- Smart resume capability (mid-step resumption using GROMACS checkpoints `.cpt`)
- Dedicated GROMACS command logging (`GMX_LOG_FILE`)
- Support for binary/ternary solvent mixtures via `solvent.conf` evaporation rates

### Version 2.0
- Modular, function-based workflow
- Automatic MDP file generation
- Independent time parameters for each simulation phase
- Comprehensive error handling
- Improved documentation
- Basic cycle skipping

### Version 1.0 (Legacy)
- Original evaporation workflow
- Template-based MDP files
- Basic automation

---

**CastingMD** - Enabling realistic polymer film formation simulations through automated solvent evaporation protocols.
