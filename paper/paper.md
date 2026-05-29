---
title: 'CastingMD: Automated Workflow for Molecular Dynamics Simulation of Polymer Film Formation via Solvent Evaporation'
tags:
  - molecular dynamics
  - polymer films
  - solvent evaporation
  - GROMACS
  - computational materials
authors:
  - name: Ricardo Stefani
    affiliation: "1, 2"
    corresponding: true
  - name: Rudiere Luz Araújo
    affiliation: 2
affiliations:
  - index: 1
    name: Department of Chemistry, Federal University of Ouro Preto, Ouro Preto, MG, Brazil
  - index: 2
    name: Graduate Program in Materials Science, Campus Araguaia, Federal University of Mato Grosso, Barra do Garças, MT, Brazil
date: 29 May 2026
bibliography: paper.bib
---

# Summary

CastingMD is a computational workflow that automates molecular dynamics (MD) simulations of polymer film formation through iterative solvent evaporation. The software orchestrates GROMACS simulations, Python-based solvent removal, and system re-equilibration cycles to realistically reproduce the casting drying process at atomistic resolution. By directly mimicking experimental solvent evaporation conditions—where removal occurs primarily along the z-axis—CastingMD captures structural reorganization, densification, and chain rearrangement that are critical to final film properties. The tool is designed for researchers in materials science, polymer chemistry, and food engineering who need to model film formation without extensive scripting expertise.

# Statement of Need

Solvent evaporation during solution casting is a critical but poorly understood determinant of polymer film structure and properties. Experimental observations show that drying rate, solvent type, and polymer-solvent interactions influence final film morphology, mechanical properties, barrier function, and bioactivity [@ahmed2017]. Molecular dynamics (MD) simulations offer unprecedented atomistic insight into this process, yet existing approaches present substantial barriers to adoption:

1. **Manual scripting burden**: Researchers typically assemble films using static packing tools (e.g., PACKMOL) or docking methods, which ignore solvent dynamics entirely. Simulating realistic evaporation requires manual coordination of repeated MD cycles, solvent removal, re-equilibration, and file management—a task prone to error and requiring significant coding expertise.

2. **Directional inaccuracy**: Prior tools (e.g., GenEvaPa) often remove solvent along all axes (x, y, z), which does not reflect real casting where evaporation is primarily unidirectional (z-axis). This fundamental mismatch compromises realism and applicability to oriented coating systems.

3. **System inflexibility**: Existing methods are typically hardcoded for specific chemical systems or thin-film geometries, making adaptation to new polymers, mixed-solvent systems, or bioactive additives laborious or infeasible.

4. **Computational reproducibility**: The absence of standardized, open-source workflows hampers reproducibility and comparison across studies.

CastingMD addresses these barriers by providing a **fully automated, adaptable, open-source framework** that removes solvent realistically along the z-axis only, supports arbitrary polymer-solvent systems including mixed-solvent formulations, manages entire simulation pipelines with minimal user intervention, and produces computationally reproducible, publication-ready workflows.

# Software Description

## Architecture and Workflow

CastingMD comprises three integrated components:

1. **evaporate.py** (Python): Parses MD trajectories, identifies solvent molecules above a z-axis threshold, removes them from both coordinate and topology files, and automatically updates system definitions.

2. **run_evaporate_loop.sh** (Bash): Master control script that orchestrates the complete workflow—it manages MD simulations, calls the evaporation script, handles re-minimization and re-equilibration, checks box dimensions for stability, and automates cycle management.

3. **evaporate.conf** (Configuration): User-editable parameters controlling simulation times (NVT, NPT, MD, annealing), temperature, pressure, MD integration parameters, number of cycles, and input file paths.

4. **solvent.conf** (Solvent Configuration): User-editable parameters controlling mixed-solvent evaporation rates.

## Workflow Cycle

Each evaporation cycle follows this sequence:

1. **Production MD**: Short MD simulation (configurable, default 1–10 ns) allowing thermal relaxation under the current solvent composition.
2. **Solvent Removal**: Python script identifies solvent molecules farthest from the polymer along the z-axis and removes them.
3. **System Re-equilibration**: Energy minimization followed by NVT (100 ps) and NPT (500 ps) equilibration to accommodate the newly solvated state.
4. **Box Adjustment**: Script checks z-box dimension and adjusts if necessary to prevent simulation instabilities.
5. **Iteration**: Steps 1–4 repeat until 99% of solvent is removed.
6. **Final Annealing**: Production run cooling the system from initial to target temperature (e.g., 308 K → 298 K over 2–10 ns) to yield a relaxed final structure.

## Key Features

- **Mixed-solvent support**: Handles systems with multiple solvents (water + acetic acid, water + ethanol, etc.).
- **Modular design**: Easy adaptation to alternative MD engines (OpenMM, AMBER) via the bash framework.
- **Automatic configuration management**: Python script creates GROMACS parameter files (.mdp) from a single configuration file, eliminating manual template editing.
- **Resume capability**: Automatically detects completed cycles and resumes interrupted simulations.
- **Comprehensive logging**: Detailed log files track all operations, enabling debugging and audit trails.

## Computational Efficiency

Typical workflow timings on a single GPU node (NVIDIA RTX 3050 equivalent):

- **System**: 4 PVA chains + 2 quercetin + ~16,500 water molecules (~20,000 atoms)
- **Per cycle**: 30–50 minutes wall-clock time (5 ns MD + equilibration)
- **Full workflow** (20 cycles): ~10–29 hours wall-clock
- **Total MD simulation time**: ~110 ns (20 cycles × 5 ns production + equilibration + annealing)

This efficiency enables iterative design exploration on modest computational resources (single GPU or HPC node with GPU acceleration).

# Validation and Results

CastingMD was validated on two polymer-solvent systems representing distinct chemical scenarios:

1. **System 1 (PVA + Quercetin)**: 4 PVA chains (30 monomers each), 2 quercetin molecules, 16,527 water molecules. Represents bioactive film with natural antioxidant.

2. **System 2 (Chitosan + Glycerol + Acetic Acid)**: 20 chitosan chains (25 monomers each), 100 glycerol molecules, 512 acetic acid molecules (~2%), 27,515 water molecules. Represents mixed-solvent system with plasticizer.

Structural metrics (RMSD, radius of gyration, density profiles, radial distribution functions) confirmed realistic densification behavior:

- **Radius of Gyration (Rg) reduction**: System 1 contracted from 4.3–4.5 nm to 1.3–1.4 nm (~70% densification), with stepwise drops synchronized to solvent removal events. System 2 showed comparable behavior (3.85–4.1 nm → 2.95–3.0 nm).

- **Density evolution**: Final density profiles showed compact, homogeneous polymer packing with elimination of solvent-rich voids, increasing from ~0.8 g/cm³ (solvated) to >1.2 g/cm³ (desiccated).

- **Molecular packing**: RDF analysis revealed transformation from solvent-separated chains to tightly packed, entangled networks with reduced long-range correlations—behavior consistent with experimental observations of cast films.

Both systems converged to stable, compact final conformations with reduced free volume—behavior expected in real casting and absent in static packing approaches.

# Comparison to Existing Software

| Feature | GenEvaPa [@harris2023] | PyThinFilm [@stroet2023] | PACKMOL [@martinez2009] | CastingMD |
|---------|----------|----------|---------|-----------|
| Automated evaporation workflow | + | + | - | + |
| Z-axis-only evaporation | - | - | N/A | + |
| Mixed-solvent support | - | + | Limited | + |
| Vacuum deposition | - | + | - | - |
| Arbitrary polymer systems | - | + | + | + |
| Open-source, active maintenance | + | + | + | + |
| No expert GROMACS required | - | + | + | + |
| Automatic re-equilibration | Partial | + | - | + |

**GenEvaPa** provides similar automation but is restricted to thin-film systems and uses multidirectional evaporation. **PyThinFilm** is a sophisticated framework supporting both vacuum deposition and solution processing with explicit solvent-vacuum interfaces and computational acceleration methods; however, it does not restrict evaporation to the z-axis and is optimized for organic semiconductors (OLEDs/OSCs). **PACKMOL** is powerful for system construction but does not simulate evaporation dynamics. **CastingMD** fills a complementary niche by emphasizing directional (z-axis only) evaporation, automatic z-box adjustment, and accessibility for solution-cast polymer-based systems (food science, biomaterials, bioactive formulations).

# Use Cases and Impact

CastingMD enables new research directions across multiple domains:

1. **Active and smart food packaging**: Researchers can now simulate how antioxidants, antimicrobial agents, or indicator compounds distribute during film drying, predicting release kinetics and effectiveness.

2. **Drug delivery systems**: Rational design of polymer-based delivery films with controlled morphology and drug localization, critical for controlled-release applications [@huo2016; @pal2006].

3. **Biodegradable membranes and coatings**: Structure-property relationships in biocompatible films for environmental and biomedical applications.

4. **Process optimization**: Quantitative prediction of how solvent type, evaporation kinetics, and additives affect final film structure—reducing costly experimental trial-and-error cycles.

The workflow is generating active research outputs and is being adopted by collaborators in materials science and polymer engineering. A full validation study is available as a preprint [@araujo2025].

# Availability and Documentation

- **Repository**: https://github.com/rstefaniufmt/castingMD
- **Software Archive**: https://doi.org/10.5281/zenodo.16728915 (Zenodo with persistent DOI)
- **Documentation**: Complete README with installation instructions, comprehensive tutorial PDF, annotated configuration templates, and working example directories with input files and expected outputs
- **License**: MIT (permissive open-source)
- **Requirements**: GROMACS (≥2020.1), Python 3.6+, MDAnalysis library, Bash shell

Installation and usage instructions are provided in the GitHub repository README. The software has been tested on Linux systems (Ubuntu 20.04+, CentOS 7+) and macOS with standard HPC cluster environments and single-GPU workstations.

# Author Contributions

R.S. conceived the software architecture, developed core algorithms, designed validation experiments, and led manuscript preparation. R.L.A. conducted validation simulations, contributed to methodology refinement, and co-authored the preprint. Both authors approve the final manuscript.

# Acknowledgments

The authors acknowledge computational resources from the Federal University of Ouro Preto. We thank the research community for feedback during development.

# References
