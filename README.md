# DFT Workshop: Quantum ESPRESSO (QE)

This repository contains the hands-on material for a Density Functional Theory (DFT) workshop based on **Quantum ESPRESSO** (PWscf), including installation notes, tutorial exercises, input templates, helper scripts, and reference outputs.

Later, this repository may be extended with material for other electronic-structure codes, but the current content is focused on Quantum ESPRESSO.

## Quick start (suggested workshop flow)

1. Choose an installation route:
   - Quantum Mobile VM (Windows): [Installation/Quantum_Mobile_Installation.md](Installation/Quantum_Mobile_Installation.md)
   - Quantum ESPRESSO on Windows via WSL: [Installation/QE_in_WSL_Windows.md](Installation/QE_in_WSL_Windows.md)
   - Remote login tools (optional): [Installation/MobaXerm_Putty_Installation.md](Installation/MobaXerm_Putty_Installation.md)
   - GUI/visualization tools (optional): [Installation/Burai_Install_Windows_Linux.md](Installation/Burai_Install_Windows_Linux.md), [Installation/VESTA_Installation.md](Installation/VESTA_Installation.md)

2. Start with the tutorials: [Tutorials_QE/](Tutorials_QE/)
   - Tutorial 1 (Linux + visualization): [Tutorials_QE/Tutorial-1/](Tutorials_QE/Tutorial-1/)
   - Tutorial 2 (SCF + convergence tests in Si): [Tutorials_QE/Tutorial-2/](Tutorials_QE/Tutorial-2/)
   - Tutorial 3 (post-processing + supercells + bands): [Tutorials_QE/Tutorial-3/](Tutorials_QE/Tutorial-3/)

3. Pick or download pseudopotentials:
   - Pseudopotentials folder + trusted download links: [psudo/README.md](psudo/README.md)

4. Verify your QE executables (examples used in this repo):
   - `pw.x` (SCF/NSCF/Bands calculations)
   - `pp.x` (charge density post-processing)
   - `dos.x` (density of states)
   - `bands.x` + `plotband.x` (band structure collection/plotting)
   - `ev.x` (equation-of-state fit used in Tutorial 2)

## Repository layout

### Installation guides

All installation guides are under [Installation/](Installation/):

- Quantum Mobile VM (recommended for a ready-to-use Linux environment): [Installation/Quantum_Mobile_Installation.md](Installation/Quantum_Mobile_Installation.md)
- QE on Windows via WSL (Ubuntu + conda): [Installation/QE_in_WSL_Windows.md](Installation/QE_in_WSL_Windows.md)
- BURAI (QE GUI): [Installation/Burai_Install_Windows_Linux.md](Installation/Burai_Install_Windows_Linux.md)
- VESTA (structure/charge visualization): [Installation/VESTA_Installation.md](Installation/VESTA_Installation.md)
- MobaXterm / PuTTY (SSH + file transfer): [Installation/MobaXerm_Putty_Installation.md](Installation/MobaXerm_Putty_Installation.md)

### Tutorials and exercises (Quantum ESPRESSO)

Top-level helper notes:
- Basic Linux commands: [Tutorials_QE/Basic_command/Basic_Linux_command.md](Tutorials_QE/Basic_command/Basic_Linux_command.md)
- Basic `vi/vim`: [Tutorials_QE/Basic_command/Basic_Vi_command.md](Tutorials_QE/Basic_command/Basic_Vi_command.md)
- Basic `awk`: [Tutorials_QE/Basic_command/Basic_awk_command.md](Tutorials_QE/Basic_command/Basic_awk_command.md)
- Basic gnuplot usage: [Tutorials_QE/Basic_command/GNU_Plot_Basic.md](Tutorials_QE/Basic_command/GNU_Plot_Basic.md)

#### Tutorial 1: Linux basics, XCrySDen, XMGrace

- Linux basics (short exercise included): [Tutorials_QE/Tutorial-1/LINUX_BASICS](Tutorials_QE/Tutorial-1/LINUX_BASICS)
- XCrySDen visualization walkthrough: [Tutorials_QE/Tutorial-1/XCRYSDEN_Tutorial](Tutorials_QE/Tutorial-1/XCRYSDEN_Tutorial)
  - Example structure inputs: [Tutorials_QE/Tutorial-1/XCRYSDEN/](Tutorials_QE/Tutorial-1/XCRYSDEN/)
- XMGrace plotting walkthrough: [Tutorials_QE/Tutorial-1/XMGRACE_Tutorial](Tutorials_QE/Tutorial-1/XMGRACE_Tutorial)
  - Example data files: [Tutorials_QE/Tutorial-1/XMGRACE/](Tutorials_QE/Tutorial-1/XMGRACE/)

#### Tutorial 2: SCF calculations and convergence tests (Si)

Main tutorial text: [Tutorials_QE/Tutorial-2/Instructions](Tutorials_QE/Tutorial-2/Instructions)

- Exercise 1: plane-wave cutoff convergence
  - Sample input: [Tutorials_QE/Tutorial-2/Exercise-1/Si.sample.in](Tutorials_QE/Tutorial-2/Exercise-1/Si.sample.in)
  - Automation script: [Tutorials_QE/Tutorial-2/Exercise-1/Si.ecutwfc.sh](Tutorials_QE/Tutorial-2/Exercise-1/Si.ecutwfc.sh)
  - Reference inputs/outputs and example data: [Tutorials_QE/Tutorial-2/Exercise-1/Reference/](Tutorials_QE/Tutorial-2/Exercise-1/Reference/)
- Exercise 2: k-point convergence
  - Automation script: [Tutorials_QE/Tutorial-2/Exercise-2/Si.sample_K.sh](Tutorials_QE/Tutorial-2/Exercise-2/Si.sample_K.sh)
  - Reference inputs/outputs and example data: [Tutorials_QE/Tutorial-2/Exercise-2/Reference/](Tutorials_QE/Tutorial-2/Exercise-2/Reference/)
- Exercise 3: lattice-constant scan and equation-of-state fit
  - Automation script: [Tutorials_QE/Tutorial-2/Exercise-3/Si.sample_alat.sh](Tutorials_QE/Tutorial-2/Exercise-3/Si.sample_alat.sh)
  - Reference inputs/outputs and example data: [Tutorials_QE/Tutorial-2/Exercise-3/Reference/](Tutorials_QE/Tutorial-2/Exercise-3/Reference/)

#### Tutorial 3: post-processing, supercells, bands

- Exercise 1: charge density extraction and DOS
  - Instructions: [Tutorials_QE/Tutorial-3/Exercise-1/INSTRUCTIONS](Tutorials_QE/Tutorial-3/Exercise-1/INSTRUCTIONS)
  - Input templates: [Tutorials_QE/Tutorial-3/Exercise-1/Si.scf.in](Tutorials_QE/Tutorial-3/Exercise-1/Si.scf.in), [Tutorials_QE/Tutorial-3/Exercise-1/Si.pp_rho.in](Tutorials_QE/Tutorial-3/Exercise-1/Si.pp_rho.in), [Tutorials_QE/Tutorial-3/Exercise-1/Si.pp_rho.2D.in](Tutorials_QE/Tutorial-3/Exercise-1/Si.pp_rho.2D.in), [Tutorials_QE/Tutorial-3/Exercise-1/dos.in](Tutorials_QE/Tutorial-3/Exercise-1/dos.in)
  - Example outputs (xsf, dos, logs): [Tutorials_QE/Tutorial-3/Exercise-1/Reference/](Tutorials_QE/Tutorial-3/Exercise-1/Reference/)
- Exercise 2: supercell calculations
  - Instructions: [Tutorials_QE/Tutorial-3/Exercise-2/INSTRUCTIONS](Tutorials_QE/Tutorial-3/Exercise-2/INSTRUCTIONS)
  - Input templates: [Tutorials_QE/Tutorial-3/Exercise-2/Si.supercell_1.in](Tutorials_QE/Tutorial-3/Exercise-2/Si.supercell_1.in), [Tutorials_QE/Tutorial-3/Exercise-2/Si.supercell_2.in](Tutorials_QE/Tutorial-3/Exercise-2/Si.supercell_2.in)
- Exercise 3: band structure (reference outputs included)
  - Instructions: [Tutorials_QE/Tutorial-3/Exercise-3/INSTRUCTIONS](Tutorials_QE/Tutorial-3/Exercise-3/INSTRUCTIONS)
  - Bands collection input: [Tutorials_QE/Tutorial-3/Exercise-3/bands.in](Tutorials_QE/Tutorial-3/Exercise-3/bands.in)
  - Brillouin-zone image: [Tutorials_QE/Tutorial-3/Exercise-3/FCC.Brillouinzone_image.png](Tutorials_QE/Tutorial-3/Exercise-3/FCC.Brillouinzone_image.png)
  - Example outputs and plotting files: [Tutorials_QE/Tutorial-3/Exercise-3/Reference/](Tutorials_QE/Tutorial-3/Exercise-3/Reference/)

## Pseudopotentials included

See [psudo/README.md](psudo/README.md) for:
- the pseudopotentials shipped with this repo
- QE-compatible pseudopotential libraries (trusted sources) and download links

## Important notes before running calculations

- Many input files and scripts contain hard-coded `pseudo_dir` paths that point to specific workshop machines. Update `pseudo_dir` to your local pseudopotential folder (for example, the included UPFs are in [psudo/](psudo/)).
- Some templates use pseudopotential filenames that do not match the UPF filenames shipped in this repo. Make sure the `ATOMIC_SPECIES` line matches an existing UPF file.
- Keep `prefix` and `outdir` consistent across SCF/NSCF and post-processing steps (`pp.x`, `dos.x`, `bands.x`), and create the `outdir` directory if the input points to a non-existent folder (for example, [Tutorials_QE/Tutorial-3/Exercise-1/dos.in](Tutorials_QE/Tutorial-3/Exercise-1/dos.in) uses `outdir = './tmp/'`).
- Each `Reference/` folder contains example inputs/outputs that can be used to validate that your run is behaving as expected.
