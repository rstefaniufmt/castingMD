# Changelog and Bug Fixes (CastingMD)

This document records all recent changes, optimizations, and bug fixes implemented in the Molecular Dynamics Film Casting simulation scripts (`evaporate_optimized.sh`, `evaporate.conf`, and auxiliary scripts).

## [2.1] - Current Version

---

## 🚀 New Features and Optimizations

### 1. Smart Resume Capability
**Script:** `evaporate_optimized.sh`
- **Description:** Implemented the ability for the script to resume the simulation from the exact point it stopped in case of an abrupt interruption (e.g., power outage or manual user cancellation).
- **How it works:** The script now checks for the presence of checkpoint (`.cpt`) and `.tpr` files. If they exist and the step has not been completed, it executes the command `gmx mdrun -cpi <file>.cpt -append`, appending the new data to the original output files instead of restarting the entire step from scratch.
- **Scope:** Applied to Energy Minimization, NVT Equilibration, NPT Equilibration, Production MD, and Final Simulation (Annealing) steps.

### 2. Improved GROMACS Logging
**Scripts:** `evaporate_optimized.sh` and `evaporate.conf`
- **Description:** The GROMACS output (which used to be suppressed or pollute the terminal) is now entirely redirected to a dedicated log file.
- **How it works:** The variable `GMX_LOG_FILE` was created (defined in `evaporate.conf` as `gromacs_output.log`). All `gmx grompp`, `gmx mdrun`, `gmx editconf`, and `gmx trjconv` commands send their standard output (stdout) and errors (stderr) to this file (`>> "$GMX_LOG_FILE" 2>&1`).
- **Benefit:** It immensely facilitates debugging in case of simulation failures, such as overlapping atoms or topology issues, keeping the terminal clean with only the script's progress (`evaporate.log`).

### 3. Support for Solvent Mixtures (Binary/Ternary)
**Script:** `evaporate.py`
- **Description:** The Python script has been updated to properly handle the evaporation of complex solvent mixtures (binary or ternary).
- **How it works:** If a `solvent.conf` file is present in the directory, the system will read the specific evaporation rates of each solvent from it.
- **Format:** The file must follow the format `SOLVENT_NAME = evaporation_rate`. For example:
  ```text
  SOL = 1.0 # (Water)
  ACTN = 1.4 # (ACETONE) 
  ```
  This allows for much more refined control over how different solvent molecules leave the system during the *casting* process.

---

## 🐛 Bug Fixes

### 1. Infinite Loop / Restarting Partially Completed Steps
**Previous Issue:** The old check in the main loop only verified if the final file of the cycle was present. If power was lost in the middle of a step (e.g., step 3 in production MD), the script wouldn't identify the partial files, rewrote the `.tpr` file by running `grompp` again, and started running MD from scratch.
**Fix:**
- Added individual checks for the presence of `.gro` files in *every* substep within `run_equilibration` and the final *annealing*.
- If the substep (e.g., NVT) already has a corresponding finished `.gro` file, the script skips it (`Skipping...`). If it doesn't, but it has a `.cpt`, it resumes.

### 2. Configuration Parser Fix
**Script:** `prepare_evaporation.py` / `evaporate.conf`
**Previous Issue:** The Python scripts could ignore custom values defined by the user in `evaporate.conf`, using "hardcoded" (native) values defined in the code instead of reading from the configuration file.
**Fix:** The script has been adjusted to correctly load and respect all constants and variables defined inside `evaporate.conf`.

---

## 📝 How to use the new features

> [!TIP]
> **Intentionally interrupting the simulation**
> You can press `Ctrl + C` in the terminal at any time to stop the script. To continue, just run `evaporate_optimized.sh` again in the same folder. It will read the files and resume from the exact point it stopped.

> [!NOTE]
> **Monitoring Errors**
> If the script reports an error and stops (e.g., "grompp failed (EM)"), open the `gromacs_output.log` file. The last lines of this file will contain the exact error returned by GROMACS.
