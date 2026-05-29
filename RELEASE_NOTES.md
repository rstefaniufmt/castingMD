# Release Notes - CastingMD v2.1

We are excited to announce the release of **CastingMD** (Molecular Dynamics Film Casting) version 2.1! This update brings significant improvements in simulation resilience, better log management, and advanced support for solvent mixtures.

Below are the details of what's new and the bug fixes implemented in this release.

---

## ✨ What's New

### 1. Smart Resume Capability
You no longer have to worry about losing your progress if a step gets interrupted by a power outage or a manual cancellation (`Ctrl + C`).
* **How it works:** The `evaporate_optimized.sh` script has been enhanced to automatically detect checkpoint (`.cpt`) and `.tpr` files. If the simulation is interrupted, running the script again will trigger an *append* (adding new data to existing files) and seamlessly resume the simulation from the exact point it stopped.
* **Where it applies:** Energy Minimization, Equilibration (NVT and NPT), Production MD, and Final Simulation (Annealing).

### 2. Improved GROMACS Logging
No more terminal clutter during GROMACS steps!
* **What changed:** All output (stdout and stderr) from GROMACS commands (`grompp`, `mdrun`, `editconf`, `trjconv`) is now entirely redirected to a dedicated log file. This keeps the terminal clean, only displaying the main script's progress (`evaporate.log`).
* **How to access it:** The log is saved in the file defined by the `GMX_LOG_FILE` variable (which defaults to `gromacs_output.log`). This greatly simplifies troubleshooting failures, topology issues, or overlapping atoms by allowing you to quickly check the exact error returned by GROMACS.

### 3. Support for Solvent Mixtures
Greater physicochemical control during the *casting* step!
* **What changed:** The `evaporate.py` script now supports the evaporation of complex solvent mixtures, such as binary or ternary systems.
* **How it works:** Simply create a `solvent.conf` file in your simulation directory and declare the specific evaporation rates for each type of molecule.
* **Usage Example (`solvent.conf`):**
  ```text
  SOL = 1.0   # Water
  ACTN = 1.4  # Acetone
  ```

---

## 🐛 Bug Fixes

### 1. Infinite Loop and Partial Step Restart Fix
* **The Issue:** The main verification loop only checked for the final file of a given cycle. If the simulation was interrupted in the middle of a step (e.g., NVT) without completing the final `.gro` file, the script would run `grompp` again and restart that MD step from scratch.
* **The Solution:** We added individual validations for *each* substep. If a substep already has its corresponding finished `.gro` file, the script simply skips it (`Skipping...`). If the `.gro` file is missing but a `.cpt` is present, the Smart Resume feature kicks in.

### 2. Configuration Parser Fix (`evaporate.conf`)
* **The Issue:** The Python scripts (`prepare_evaporation.py`) were occasionally ignoring custom values set by the user in `evaporate.conf`, defaulting back to hardcoded values defined in the source code.
* **The Solution:** The parser was adjusted to ensure that all variables and constants defined in `evaporate.conf` are properly prioritized, loaded, and respected during execution.

---

> **Pro Tip:** To test the new resume functionality, feel free to pause your simulation with `Ctrl + C` at any time and run the main script again. CastingMD will pick up right where it left off without losing any data!
