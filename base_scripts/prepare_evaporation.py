#!/usr/bin/env python3

"""
PREPARE EVAPORATION - MDP File Generator
Generates GROMACS MDP files for film casting MD simulations
Reads parameters from evaporate.conf and creates all necessary input files
"""

import os
import sys
from pathlib import Path

# ============================================================================
# CONFIGURATION
# ============================================================================

CONFIG_FILE = "evaporate.conf"
DEFAULT_PARAMS = {
    # Simulation times
    'NVT_PS': 100,         # NVT equilibration time (ps)
    'NPT_PS': 1000,        # NPT equilibration time (ps)
    'SEG_NS': 1,           # Production MD time per cycle (ns)
    'ANNEAL_NS': 10,       # Annealing time (ns)
    
    # Thermodynamic parameters
    'TEMP': 308,           # Temperature (K)
    'PRESSURE': 1.0,       # Pressure (bar)
    'ANNEAL_TEMP': 298,    # Final annealing temperature (K)
    
    # MD parameters
    'DT': 0.001,           # Time step (ps)
    'RLIST': 1.2,          # Neighbor list cutoff (nm)
    'RVDW': 1.2,           # VdW cutoff (nm)
    'RCOULOMB': 1.2,       # Coulomb cutoff (nm)
    
    # Output frequencies (in steps, relative to dt)
    'NVT_OUT_FREQ': 5000,  # NVT output every 5000 steps
    'NPT_OUT_FREQ': 500,   # NPT output every 500 steps
    'MD_OUT_FREQ': 100000, # MD output every 100000 steps (100 ps for dt=0.001)
}

# ============================================================================
# FUNCTIONS
# ============================================================================

def load_config(config_file):
    """
    Load parameters from evaporate.conf
    Returns dictionary with configuration parameters
    """
    params = DEFAULT_PARAMS.copy()
    
    if not os.path.exists(config_file):
        print(f"[WARNING] Configuration file '{config_file}' not found.")
        print(f"[INFO] Using default parameters: {params}")
        return params
    
    print(f"[INFO] Loading configuration from '{config_file}'...")
    
    with open(config_file, 'r') as f:
        for line in f:
            line = line.strip()
            # Skip comments and empty lines
            if not line or line.startswith('#'):
                continue
            
            # Parse variable assignments (format: VAR=value or VAR="value")
            if '=' in line:
                # Remove 'export' if present
                line = line.replace('export', '').strip()
                
                key, value = line.split('=', 1)
                key = key.strip()
                value = value.strip().strip('"').strip("'")
                
                # Convert known numeric parameters
                numeric_params = [
                    'SEG_NS', 'NVT_PS', 'NPT_PS', 'ANNEAL_NS',
                    'TEMP', 'PRESSURE', 'ANNEAL_TEMP',
                    'DT', 'RLIST', 'RVDW', 'RCOULOMB',
                    'NVT_OUT_FREQ', 'NPT_OUT_FREQ', 'MD_OUT_FREQ'
                ]
                
                if key in numeric_params:
                    try:
                        params[key] = float(value)
                    except ValueError:
                        print(f"[WARNING] Could not convert {key}={value} to number. Using default.")
    
    print(f"[INFO] Loaded parameters: {params}")
    return params


def generate_em_mdp():
    """Generate energy minimization MDP file"""
    
    content = """; LINES STARTING WITH ';' ARE COMMENTS
title		    = Minimization	; Title of run

; Parameters describing what to do, when to stop and what to save
integrator	    = steep		; Algorithm (steep = steepest descent minimization)
emtol		    = 500.0  	; Stop minimization when the maximum force < 1.0 kJ/mol
emstep          = 0.01      ; Energy step size
nsteps		    = 10000	  	; Maximum number of (minimization) steps to perform

; Parameters describing how to find the neighbors of each atom and how to calculate the interactions
constraints     = none
nstlist		    = 1		        ; Frequency to update the neighbor list and long range forces
cutoff-scheme   = Verlet
ns_type		    = grid		    ; Method to determine neighbor list (simple, grid)
rlist		    = 1.2		    ; Cut-off for making neighbor list (short range forces)
coulombtype	    = PME		    ; Treatment of long range electrostatic interactions
rcoulomb	    = 1.2		    ; long range electrostatic cut-off
vdwtype         = cutoff
vdw-modifier    = force-switch
rvdw-switch     = 1.0
rvdw		    = 1.2		    ; long range Van der Waals cut-off
pbc             = xyz 		    ; Periodic Boundary Conditions
DispCorr        = no
"""
    
    return content


def generate_nvt_mdp(params):
    """Generate NVT equilibration MDP file"""
    
    nvt_ps = params['NVT_PS']
    temp = params['TEMP']
    dt = params['DT']
    rlist = params['RLIST']
    rvdw = params['RVDW']
    rcoulomb = params['RCOULOMB']
    out_freq = int(params['NVT_OUT_FREQ'])
    
    # Calculate nsteps based on NVT_PS
    nsteps = int(nvt_ps / dt)
    
    content = f"""; NVT Equilibration
; Generated by prepare_evaporation.py

title                   = NVT Equilibration
define                  = -DPOSRES      ; Position restrain heavy atoms

; Run parameters
integrator              = md            ; Leap-frog integrator
nsteps                  = {nsteps}      ; {nvt_ps:.1f} ps
dt                      = {dt}          ; {dt*1000:.1f} fs

; Output control
nstenergy               = {out_freq}    ; Save energies every {out_freq} steps
nstlog                  = {out_freq}    ; Update log every {out_freq} steps
nstxout-compressed      = {out_freq}    ; Save coordinates every {out_freq} steps

; Bond parameters
continuation            = yes           ; Continuation from EM
constraint_algorithm    = lincs         ; Holonomic constraints
constraints             = h-bonds       ; Constrain H-bonds
lincs_iter              = 2             ; Accuracy of LINCS
lincs_order             = 4             ; Also related to accuracy

; Neighbor searching and VdW
cutoff-scheme           = Verlet
ns_type                 = grid          ; Search neighboring grid cells
nstlist                 = 20            ; Largely irrelevant with Verlet
rlist                   = {rlist}
vdwtype                 = cutoff
vdw-modifier            = force-switch
rvdw-switch             = 1.0
rvdw                    = {rvdw}        ; VdW cutoff (nm)

; Electrostatics
coulombtype             = PME           ; Particle Mesh Ewald
rcoulomb                = {rcoulomb}    ; Coulomb cutoff (nm)
pme_order               = 4             ; Cubic interpolation
fourierspacing          = 0.16          ; Grid spacing for FFT

; Temperature coupling
tcoupl                  = V-rescale     ; Modified Berendsen thermostat
tc-grps                 = System
tau_t                   = 1.0           ; Time constant (ps)
ref_t                   = {temp}        ; Reference temperature (K)

; Pressure coupling
pcoupl                  = no            ; No pressure coupling in NVT

; Periodic boundary conditions
pbc                     = xyz           ; 3-D PBC

; Dispersion correction
DispCorr                = no

; Velocity generation
gen_vel                 = no            ; Velocities from previous step
gen_temp                = {temp}        ; Temperature for Maxwell distribution
gen_seed                = -1            ; Random seed
"""
    
    return content


def generate_npt_mdp(params):
    """Generate NPT equilibration MDP file"""
    
    npt_ps = params['NPT_PS']
    temp = params['TEMP']
    pressure = params['PRESSURE']
    dt = params['DT']
    rlist = params['RLIST']
    rvdw = params['RVDW']
    rcoulomb = params['RCOULOMB']
    out_freq = int(params['NPT_OUT_FREQ'])
    
    # Calculate nsteps based on NPT_PS
    nsteps = int(npt_ps / dt)
    
    content = f"""; NPT Equilibration
; Generated by prepare_evaporation.py

title                   = NPT Equilibration
define                  = -DPOSRES      ; Position restrain heavy atoms

; Run parameters
integrator              = md            ; Leap-frog integrator
nsteps                  = {nsteps}      ; {npt_ps:.1f} ps ({npt_ps/1000:.1f} ns)
dt                      = {dt}          ; {dt*1000:.1f} fs

; Output control
nstenergy               = {out_freq}    ; Save energies every {out_freq} steps
nstlog                  = {out_freq}    ; Update log every {out_freq} steps
nstxout-compressed      = {out_freq}    ; Save coordinates every {out_freq} steps

; Bond parameters
continuation            = yes           ; Continuing from NVT
constraint_algorithm    = lincs         ; Holonomic constraints
constraints             = h-bonds       ; Constrain H-bonds
lincs_iter              = 1             ; Accuracy of LINCS
lincs_order             = 4             ; Also related to accuracy

; Neighbor searching and VdW
cutoff-scheme           = Verlet
ns_type                 = grid          ; Search neighboring grid cells
nstlist                 = 20            ; Largely irrelevant with Verlet
rlist                   = {rlist}
vdwtype                 = cutoff
vdw-modifier            = force-switch
rvdw-switch             = 1.0
rvdw                    = {rvdw}        ; VdW cutoff (nm)

; Electrostatics
coulombtype             = PME           ; Particle Mesh Ewald
rcoulomb                = {rcoulomb}    ; Coulomb cutoff (nm)
pme_order               = 4             ; Cubic interpolation
fourierspacing          = 0.16          ; Grid spacing for FFT

; Temperature coupling
tcoupl                  = V-rescale     ; Modified Berendsen thermostat
tc-grps                 = System
tau_t                   = 1.0           ; Time constant (ps)
ref_t                   = {temp}        ; Reference temperature (K)

; Pressure coupling
pcoupl                  = C-rescale     ; Pressure coupling on for NPT
pcoupltype              = semiisotropic ; Semiisotropic scaling
tau_p                   = 2.0           ; Time constant (ps)
ref_p                   = 0.0 {pressure} ; Reference pressure (bar)
compressibility         = 0.0 4.5e-5    ; Isothermal compressibility
refcoord_scaling        = com

; Periodic boundary conditions
pbc                     = xyz           ; 3-D PBC

; Dispersion correction
DispCorr                = no

; Velocity generation
gen_vel                 = no            ; Continue from NVT
"""
    
    return content


def generate_md_mdp(params):
    """Generate production MD MDP file"""
    
    seg_ns = params['SEG_NS']
    temp = params['TEMP']
    pressure = params['PRESSURE']
    dt = params['DT']
    rlist = params['RLIST']
    rvdw = params['RVDW']
    rcoulomb = params['RCOULOMB']
    out_freq = int(params['MD_OUT_FREQ'])
    
    # Calculate nsteps based on SEG_NS
    nsteps = int(seg_ns * 1000000)  # SEG_NS (ns) * 1,000,000
    
    content = f"""; Production MD Simulation
; Generated by prepare_evaporation.py

title                   = Production MD Simulation

; Run parameters
integrator              = md            ; Leap-frog integrator
nsteps                  = {nsteps}      ; {seg_ns} ns
dt                      = {dt}          ; {dt*1000:.1f} fs

; Output control
nstxout                 = {out_freq}    ; Save coordinates to .trr every {out_freq} steps
nstvout                 = {out_freq}    ; Save velocities to .trr every {out_freq} steps
nstenergy               = {out_freq}    ; Save energies every {out_freq} steps
nstlog                  = {out_freq}    ; Update log every {out_freq} steps
nstxout-compressed      = {out_freq}    ; Save compressed coords every {out_freq} steps

; Bond parameters
continuation            = yes           ; Continuing from NPT
constraint_algorithm    = lincs         ; Holonomic constraints
constraints             = h-bonds       ; Constrain H-bonds
lincs_iter              = 2             ; Accuracy of LINCS
lincs_order             = 4             ; Also related to accuracy

; Neighbor searching and VdW
cutoff-scheme           = Verlet
ns_type                 = grid          ; Search neighboring grid cells
nstlist                 = 20            ; Largely irrelevant with Verlet
rlist                   = {rlist}
vdwtype                 = cutoff
vdw-modifier            = force-switch
rvdw-switch             = 1.0
rvdw                    = {rvdw}        ; VdW cutoff (nm)

; Electrostatics
coulombtype             = PME           ; Particle Mesh Ewald
rcoulomb                = {rcoulomb}    ; Coulomb cutoff (nm)
pme_order               = 4             ; Cubic interpolation
fourierspacing          = 0.16          ; Grid spacing for FFT

; Temperature coupling
tcoupl                  = V-rescale     ; Modified Berendsen thermostat
tc-grps                 = System
tau_t                   = 1.0           ; Time constant (ps)
ref_t                   = {temp}        ; Reference temperature (K)

; Pressure coupling
pcoupl                  = C-rescale     ; Pressure coupling on
pcoupltype              = semiisotropic ; Semiisotropic scaling
tau_p                   = 5.0           ; Time constant (ps)
ref_p                   = 0.0 {pressure} ; Reference pressure (bar)
compressibility         = 0.0 4.5e-5    ; Isothermal compressibility

; Periodic boundary conditions
pbc                     = xyz           ; 3-D PBC

; Dispersion correction
DispCorr                = no

; Velocity generation
gen_vel                 = no            ; Continue from NPT
"""
    
    return content


def generate_annealing_mdp(params):
    """Generate annealing (cooling) MD MDP file"""
    
    anneal_ns = params['ANNEAL_NS']
    temp = params['TEMP']
    anneal_temp = params['ANNEAL_TEMP']
    pressure = params['PRESSURE']
    dt = params['DT']
    rlist = params['RLIST']
    rvdw = params['RVDW']
    rcoulomb = params['RCOULOMB']
    
    # Calculate nsteps based on ANNEAL_NS
    nsteps = int(anneal_ns * 1000000)  # ANNEAL_NS (ns) * 1,000,000
    anneal_time_ps = anneal_ns * 1000  # Convert ns to ps
    
    # Output frequency: every 50 ps for annealing
    out_freq = int(50 / dt)
    out_freq_energy = int(5 / dt)
    
    content = f"""; Annealing (Cooling) MD Simulation
; Generated by prepare_evaporation.py

title                   = Annealing - Cooling to {anneal_temp} K

; Run parameters
integrator              = md            ; Leap-frog integrator
nsteps                  = {nsteps}      ; {anneal_ns} ns
dt                      = {dt}          ; {dt*1000:.1f} fs

; Output control
nstxout                 = {out_freq}    ; Save coordinates every 50 ps
nstvout                 = {out_freq}    ; Save velocities every 50 ps
nstenergy               = {out_freq_energy} ; Save energies every 5 ps
nstlog                  = {out_freq_energy} ; Update log every 5 ps
nstxout-compressed      = {out_freq_energy} ; Save compressed coords every 5 ps

; Bond parameters
continuation            = yes           ; Continuing from production
constraint_algorithm    = lincs         ; Holonomic constraints
constraints             = h-bonds       ; Constrain H-bonds
lincs_iter              = 2             ; Accuracy of LINCS
lincs_order             = 4             ; Also related to accuracy

; Neighbor searching and VdW
cutoff-scheme           = Verlet
ns_type                 = grid          ; Search neighboring grid cells
nstlist                 = 20            ; Largely irrelevant with Verlet
rlist                   = {rlist}
vdwtype                 = cutoff
vdw-modifier            = force-switch
rvdw-switch             = 1.0
rvdw                    = {rvdw}        ; VdW cutoff (nm)

; Electrostatics
coulombtype             = PME           ; Particle Mesh Ewald
rcoulomb                = {rcoulomb}    ; Coulomb cutoff (nm)
pme_order               = 4             ; Cubic interpolation
fourierspacing          = 0.16          ; Grid spacing for FFT

; Temperature coupling with annealing
tcoupl                  = V-rescale     ; Modified Berendsen thermostat
tc-grps                 = System
tau_t                   = 1.0           ; Time constant (ps)
ref_t                   = {temp}        ; Will be overridden by annealing

; Simulated annealing
annealing               = single        ; Single sequence
annealing-npoints       = 2             ; Two points: start and end
annealing-time          = 0 {anneal_time_ps:.0f} ; 0 ps and {anneal_time_ps:.0f} ps ({anneal_ns} ns)
annealing-temp          = {temp} {anneal_temp} ; Cool from {temp} to {anneal_temp} K

; Pressure coupling
pcoupl                  = C-rescale     ; Pressure coupling on
pcoupltype              = semiisotropic ; Semiisotropic scaling
tau_p                   = 5.0           ; Time constant (ps)
ref_p                   = 0.0 {pressure} ; Reference pressure (bar)
compressibility         = 0.0 4.5e-5    ; Isothermal compressibility

; Periodic boundary conditions
pbc                     = xyz           ; 3-D PBC

; Dispersion correction
DispCorr                = no

; Velocity generation
gen_vel                 = no            ; Continue from production
"""
    
    return content


def write_mdp_file(filename, content):
    """Write MDP content to file"""
    with open(filename, 'w') as f:
        f.write(content)
    print(f"[SUCCESS] Created: {filename}")


def create_example_config():
    """Create an example evaporate.conf file"""
    
    config_content = """# ============================================================================
# EVAPORATION SIMULATION CONFIGURATION
# Configuration file for film casting MD simulations
# ============================================================================

# Simulation time parameters
NVT_PS=100                  # NVT equilibration time (ps)
NPT_PS=1000                 # NPT equilibration time (ps) 
SEG_NS=1                    # Production MD time per cycle (ns)
ANNEAL_NS=10                # Final annealing time (ns)
CYCLES=20                   # Total number of evaporation cycles

# Thermodynamic parameters
TEMP=308                    # Simulation temperature (K)
PRESSURE=1.0                # Simulation pressure (bar)
ANNEAL_TEMP=298             # Final annealing temperature (K)

# MD integration parameters
DT=0.001                    # Time step (ps) - 0.001 = 1 fs, 0.002 = 2 fs
RLIST=1.2                   # Neighbor list cutoff (nm)
RVDW=1.2                    # Van der Waals cutoff (nm)
RCOULOMB=1.2                # Coulomb cutoff (nm)

# Output frequencies (in steps)
NVT_OUT_FREQ=5000           # NVT output every 5000 steps (~5 ps for dt=0.001)
NPT_OUT_FREQ=500            # NPT output every 500 steps (~0.5 ps for dt=0.001)
MD_OUT_FREQ=100000          # MD output every 100000 steps (~100 ps for dt=0.001)

# Input/Output file names
TOP=film.top                # Initial topology file
GRO=npt.gro                 # Initial structure file
CPT=npt.cpt                 # Initial checkpoint file
LOG_FILE=evaporate.log      # Log file name
OUT_FILM=final_film         # Final film output name

# Software
PYTHON=python3              # Python executable with MDAnalysis
"""
    
    with open("evaporate.conf.example", 'w') as f:
        f.write(config_content)
    
    print("[INFO] Created example configuration file: evaporate.conf.example")



# ============================================================================
# MAIN
# ============================================================================

def main():
    """Main function"""
    
    print("=" * 70)
    print("PREPARE EVAPORATION - MDP File Generator")
    print("=" * 70)
    print()
    
    # Load configuration
    params = load_config(CONFIG_FILE)
    print()
    
    # Generate MDP files
    print("[INFO] Generating MDP files...")
    print()
    
    mdp_files = {
        'em.mdp': generate_em_mdp(),
        'nvt_eq.mdp': generate_nvt_mdp(params),
        'npt_eq.mdp': generate_npt_mdp(params),
        'md.mdp': generate_md_mdp(params),
        'md_anneling.mdp': generate_annealing_mdp(params)
    }
    
    # Write all files
    for filename, content in mdp_files.items():
        write_mdp_file(filename, content)
    
    print()
    print("=" * 70)
    print("MDP FILE GENERATION COMPLETED")
    print("=" * 70)
    print()
    print("Generated files:")
    for filename in mdp_files.keys():
        print(f"  ✓ {filename}")
    print()
    print(f"Key parameters:")
    print(f"  - NVT equilibration: {params['NVT_PS']} ps")
    print(f"  - NPT equilibration: {params['NPT_PS']} ps ({params['NPT_PS']/1000:.1f} ns)")
    print(f"  - Production MD per cycle: {params['SEG_NS']} ns ({int(params['SEG_NS']*1000000)} steps)")
    print(f"  - Annealing: {params['ANNEAL_NS']} ns ({params['TEMP']} K → {params['ANNEAL_TEMP']} K)")
    print(f"  - Temperature: {params['TEMP']} K")
    print(f"  - Pressure: {params['PRESSURE']} bar")
    print(f"  - Time step: {params['DT']} ps ({params['DT']*1000} fs)")
    print()
    
    # Create example config if doesn't exist
    if not os.path.exists(CONFIG_FILE):
        print(f"[INFO] Configuration file '{CONFIG_FILE}' not found.")
        create_example_config()
        print(f"[INFO] Please edit 'evaporate.conf.example' and save as '{CONFIG_FILE}'")
    
    print()
    print("Ready to run evaporation simulation!")
    print()


if __name__ == "__main__":
    main()
