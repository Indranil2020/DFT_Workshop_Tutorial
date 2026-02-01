# Quantum ESPRESSO Quick Reference Card
## Comprehensive DFT Workflow Guide

---

## 1. Complete DFT Workflow Checklist

```
[ ] Phase 1: Database Search
    [ ] Search Materials Project / OQMD / AFLOW
    [ ] Download structure file (CIF/POSCAR)
    [ ] Check literature for experimental data

[ ] Phase 2: Structure Validation
    [ ] Verify charge neutrality
    [ ] Check bond lengths are reasonable
    [ ] Visualize structure (VESTA)
    [ ] Identify space group

[ ] Phase 3: DFT Setup
    [ ] Choose XC functional (PBE default)
    [ ] Select pseudopotentials (SSSP)
    [ ] Convergence test: ecutwfc
    [ ] Convergence test: k-points

[ ] Phase 4: Ground State
    [ ] Check for magnetic elements
    [ ] Test magnetic configurations (if needed)
    [ ] Structure optimization (vc-relax or EOS)

[ ] Phase 5: Stability Tests
    [ ] Formation energy < 0
    [ ] Energy above hull < 25 meV/atom
    [ ] No imaginary phonons
    [ ] Born stability criteria satisfied

[ ] Phase 6: Properties (ONLY after stability!)
    [ ] Band structure
    [ ] DOS / PDOS
    [ ] Other properties as needed
```

---

## 2. Unit Conversions

| From       | To         | Multiply by   |
|------------|------------|---------------|
| Bohr       | Angstrom   | 0.529177      |
| Angstrom   | Bohr       | 1.889726      |
| Ry         | eV         | 13.6057       |
| eV         | Ry         | 0.073499      |
| Ry/Bohr^3  | GPa        | 14710.5       |
| GPa        | Ry/Bohr^3  | 6.798e-5      |
| Ry         | kJ/mol     | 1312.75       |
| Ry         | kcal/mol   | 313.755       |
| Ha         | eV         | 27.2114       |
| Ha         | Ry         | 2.0           |
| eV/atom    | meV/atom   | 1000          |

### Quick Mental Math
- 1 Ry ~ 13.6 eV ~ 1313 kJ/mol
- 1 Bohr ~ 0.53 Angstrom
- 1 GPa ~ 10 kbar

---

## 3. Key QE Parameters

| Parameter   | Description              | Typical Values         |
|-------------|--------------------------|------------------------|
| ecutwfc     | Wavefunction cutoff      | 40-80 Ry               |
| ecutrho     | Charge density cutoff    | 4-12x ecutwfc          |
| conv_thr    | SCF convergence          | 1e-8 Ry                |
| K_POINTS    | BZ sampling              | 6x6x6 to 12x12x12      |
| degauss     | Smearing width           | 0.01-0.02 Ry           |
| mixing_beta | Charge mixing            | 0.3-0.7                |
| forc_conv_thr | Force convergence      | 1e-4 Ry/Bohr           |
| press_conv_thr | Pressure convergence  | 0.5 kbar               |

---

## 4. Ecutrho Guidelines

| PP Type                  | ecutrho        | Notes                        |
|--------------------------|----------------|------------------------------|
| NC (Norm-conserving)     | 4x ecutwfc     | Hardest, most accurate       |
| US (Ultrasoft)           | 8-12x ecutwfc  | Softer, requires more rho    |
| PAW                      | 8x ecutwfc     | All-electron reconstruction  |

### SSSP Recommendations
- **SSSP Efficiency**: Lower cutoffs, faster calculations
- **SSSP Precision**: Higher cutoffs, more accurate

---

## 5. Calculation Types

| Type     | Purpose                    | Key Parameters                    |
|----------|----------------------------|-----------------------------------|
| scf      | Self-consistent field      | conv_thr                          |
| relax    | Ionic optimization         | forc_conv_thr, nstep              |
| vc-relax | Cell + ionic optimization  | press_conv_thr, cell_dofree       |
| bands    | Band structure             | nbnd, k-path (after scf)          |
| nscf     | Dense mesh for DOS         | nbnd, dense k-mesh (after scf)    |

### Workflow Dependencies
```
scf ─────────────────────> bands ──> bands.x ──> Plot
  │
  └───> nscf ──> dos.x ──────────────────────> DOS plot
          │
          └───> projwfc.x ─────────────────> PDOS plot
```

---

## 6. High-Symmetry Points

### FCC (Diamond, Zinc-blende, Rocksalt)

| Point | Coordinates (crystal)     |
|-------|---------------------------|
| G     | (0, 0, 0)                 |
| X     | (0.5, 0, 0.5)             |
| W     | (0.5, 0.25, 0.75)         |
| K     | (0.375, 0.375, 0.75)      |
| L     | (0.5, 0.5, 0.5)           |
| U     | (0.625, 0.25, 0.625)      |

**Common Path:** G - X - W - K - G - L - U - W - L - K | U - X

### BCC (Body-Centered Cubic)

| Point | Coordinates (crystal)     |
|-------|---------------------------|
| G     | (0, 0, 0)                 |
| H     | (0.5, -0.5, 0.5)          |
| N     | (0, 0, 0.5)               |
| P     | (0.25, 0.25, 0.25)        |

**Common Path:** G - H - N - G - P - H | P - N

### Hexagonal (HCP, Wurtzite)

| Point | Coordinates (crystal)     |
|-------|---------------------------|
| G     | (0, 0, 0)                 |
| M     | (0.5, 0, 0)               |
| K     | (1/3, 1/3, 0)             |
| A     | (0, 0, 0.5)               |
| L     | (0.5, 0, 0.5)             |
| H     | (1/3, 1/3, 0.5)           |

**Common Path:** G - M - K - G - A - L - H - A | L - M | K - H

### Simple Cubic

| Point | Coordinates (crystal)     |
|-------|---------------------------|
| G     | (0, 0, 0)                 |
| X     | (0.5, 0, 0)               |
| M     | (0.5, 0.5, 0)             |
| R     | (0.5, 0.5, 0.5)           |

**Common Path:** G - X - M - G - R - X | M - R

---

## 7. Born Stability Criteria

### Cubic Systems (3 independent constants: C11, C12, C44)

```
C11 > 0
C11 - C12 > 0  (tetragonal shear)
C11 + 2*C12 > 0  (bulk modulus > 0)
C44 > 0
```

**Derived Properties:**
- Bulk modulus: B = (C11 + 2*C12) / 3
- Shear modulus (Voigt): G_V = (C11 - C12 + 3*C44) / 5
- Shear modulus (Reuss): G_R = 5*(C11 - C12)*C44 / [4*C44 + 3*(C11 - C12)]

### Hexagonal Systems (5 independent constants: C11, C12, C13, C33, C44)

```
C11 > |C12|
C33*(C11 + C12) > 2*C13^2
C44 > 0
C66 = (C11 - C12)/2 > 0
```

### Tetragonal Systems (6-7 independent constants)

```
C11 > 0, C33 > 0, C44 > 0, C66 > 0
C11 - C12 > 0
C11 + C33 - 2*C13 > 0
2*C11 + C33 + 2*C12 + 4*C13 > 0
```

---

## 8. Shannon Ionic Radii (Selected)

### Alkali Metals (CN = 6)

| Ion  | Radius (Angstrom) |
|------|-------------------|
| Li+  | 0.76              |
| Na+  | 1.02              |
| K+   | 1.38              |
| Rb+  | 1.52              |
| Cs+  | 1.67              |

### Alkaline Earth Metals (CN = 6)

| Ion  | Radius (Angstrom) |
|------|-------------------|
| Mg2+ | 0.72              |
| Ca2+ | 1.00              |
| Sr2+ | 1.18              |
| Ba2+ | 1.35              |

### Transition Metals (CN = 6)

| Ion  | Radius (Angstrom) | Note       |
|------|-------------------|------------|
| Ti4+ | 0.605             |            |
| V5+  | 0.54              |            |
| Cr3+ | 0.615             |            |
| Mn2+ | 0.83 / 0.67       | HS / LS    |
| Fe2+ | 0.78 / 0.61       | HS / LS    |
| Fe3+ | 0.645 / 0.55      | HS / LS    |
| Co2+ | 0.745 / 0.65      | HS / LS    |
| Ni2+ | 0.69              |            |
| Cu2+ | 0.73              |            |
| Zn2+ | 0.74              |            |

### Anions (CN = 6)

| Ion  | Radius (Angstrom) |
|------|-------------------|
| O2-  | 1.40              |
| S2-  | 1.84              |
| F-   | 1.33              |
| Cl-  | 1.81              |
| Br-  | 1.96              |
| I-   | 2.20              |

**Note:** HS = High Spin, LS = Low Spin

---

## 9. Birch-Murnaghan Equation of State

### Third-Order Birch-Murnaghan EOS

```
E(V) = E0 + (9*V0*B0/16) * [(eta-1)^3 * B0' + (eta-1)^2 * (6 - 4*eta)]

where:
  eta = (V0/V)^(2/3)
  E0  = equilibrium energy
  V0  = equilibrium volume
  B0  = bulk modulus at V0
  B0' = pressure derivative of bulk modulus (typically 3.5-4.5)
```

### Python Implementation

```python
import numpy as np
from scipy.optimize import curve_fit

def birch_murnaghan(V, E0, V0, B0, B0_prime):
    """Third-order Birch-Murnaghan EOS."""
    eta = (V0 / V) ** (2.0 / 3.0)
    E = E0 + (9.0 * V0 * B0 / 16.0) * (
        (eta - 1.0) ** 3 * B0_prime +
        (eta - 1.0) ** 2 * (6.0 - 4.0 * eta)
    )
    return E

# Fit to data
# popt, pcov = curve_fit(birch_murnaghan, volumes, energies,
#                        p0=[E_guess, V_guess, B_guess, 4.0])
```

### Volume-Lattice Relations

| Structure | Primitive Volume | Conventional Volume |
|-----------|------------------|---------------------|
| FCC       | a^3 / 4          | a^3                 |
| BCC       | a^3 / 2          | a^3                 |
| SC        | a^3              | a^3                 |
| HCP       | sqrt(3)*a^2*c/2  | sqrt(3)*a^2*c       |
| Diamond   | a^3 / 4          | a^3                 |

---

## 10. Common QE Commands

### Basic Calculations

```bash
# SCF calculation
pw.x -in scf.in > scf.out

# Relaxation
pw.x -in relax.in > relax.out

# Variable-cell relaxation
pw.x -in vc-relax.in > vc-relax.out

# Band structure (after SCF)
pw.x -in bands.in > bands.out
bands.x -in bands_pp.in > bands_pp.out

# DOS (after SCF)
pw.x -in nscf.in > nscf.out
dos.x -in dos.in > dos.out

# Projected DOS
projwfc.x -in projwfc.in > projwfc.out
```

### Post-Processing

```bash
# Charge density visualization
pp.x -in pp.in > pp.out

# Phonon calculation
ph.x -in ph.in > ph.out

# Phonon post-processing
q2r.x -in q2r.in > q2r.out
matdyn.x -in matdyn.in > matdyn.out
```

### Parallel Execution

```bash
# MPI parallel
mpirun -np 4 pw.x -in scf.in > scf.out

# With threading
export OMP_NUM_THREADS=2
mpirun -np 4 pw.x -in scf.in > scf.out

# Pool parallelization (for k-points)
mpirun -np 8 pw.x -npool 4 -in scf.in > scf.out
```

---

## 11. Troubleshooting Quick Tips

| Problem                    | Likely Cause                  | Solution                              |
|----------------------------|-------------------------------|---------------------------------------|
| SCF not converging         | Mixing too aggressive         | Reduce mixing_beta to 0.2-0.3         |
| SCF oscillating            | Metallic system               | Add smearing (occupations='smearing') |
| Negative phonon frequencies| Unrelaxed structure           | Fully relax structure first           |
| Wrong band gap             | DFT limitation                | Expected; use HSE/GW for accuracy     |
| Memory error               | System too large              | Reduce ecutwfc or k-points            |
| Stress too high            | Lattice not optimized         | Run vc-relax or EOS fitting           |
| Forces not converging      | Poor initial structure        | Check atomic positions, use BFGS      |
| Charge sloshing            | Large cell / bad mixing       | Use 'local-TF' mixing_mode            |
| Slow convergence           | Complex electronic structure  | Try different diagonalization         |
| Segmentation fault         | Memory overflow               | Increase stack size: ulimit -s unlimited |

### SCF Convergence Tricks

```
# For difficult cases, add to &ELECTRONS:
  mixing_mode = 'local-TF'
  mixing_beta = 0.2
  electron_maxstep = 200
  diagonalization = 'david'
```

### Memory Estimation

```
Memory ~ (ecutwfc)^1.5 * (# k-points) * (# atoms)
```

---

## 12. Essential URLs

### Pseudopotential Libraries
- **SSSP (Standard Solid-State Pseudopotentials)**
  https://www.materialscloud.org/discover/sssp

- **PSlibrary**
  https://dalcorso.github.io/pslibrary/

- **Pseudo-Dojo**
  http://www.pseudo-dojo.org/

### Structure Databases
- **Materials Project**
  https://next-gen.materialsproject.org/

- **OQMD (Open Quantum Materials Database)**
  https://oqmd.org/

- **AFLOW**
  https://aflowlib.org/

- **ICSD (Inorganic Crystal Structure Database)**
  https://icsd.products.fiz-karlsruhe.de/

- **COD (Crystallography Open Database)**
  https://www.crystallography.net/

### Tools
- **QE Input Generator**
  https://www.materialscloud.org/work/tools/qeinputgenerator

- **SeeK-path (k-path finder)**
  https://www.materialscloud.org/work/tools/seekpath

- **VESTA (Visualization)**
  https://jp-minerals.org/vesta/

- **XCrySDen**
  http://www.xcrysden.org/

### Documentation
- **QE Official Site**
  https://www.quantum-espresso.org/

- **QE Tutorials**
  https://www.quantum-espresso.org/tutorials/

- **QE User Guide**
  https://www.quantum-espresso.org/Doc/INPUT_PW.html

---

## 13. Reference Values for Validation

### Silicon (PBE)

| Property            | DFT (PBE)     | Experiment    |
|---------------------|---------------|---------------|
| Lattice (Angstrom)  | ~5.47         | 5.43          |
| Bulk modulus (GPa)  | ~90-95        | 98.8          |
| Band gap (eV)       | ~0.5-0.6      | 1.17 (indirect)|

### Common DFT Errors

| Property        | Typical DFT Error          |
|-----------------|----------------------------|
| Lattice constant| +1-2% (PBE), -1% (LDA)     |
| Bulk modulus    | +/- 10%                    |
| Band gap        | ~40-50% underestimation    |
| Cohesive energy | ~10-20% overestimation     |

---

## 14. Convergence Testing Protocol

### ecutwfc Convergence

```
1. Start with recommended value from PP file
2. Test: ecutwfc = 30, 40, 50, 60, 70, 80 Ry
3. Converged when Delta_E < 1 meV/atom
4. Use value where convergence is achieved + 10 Ry safety margin
```

### k-point Convergence

```
1. Start with 4x4x4 mesh
2. Increase systematically: 6x6x6, 8x8x8, 10x10x10, 12x12x12
3. Converged when Delta_E < 1 meV/atom
4. For metals: use denser mesh + smearing
```

### Convergence Thresholds

| Calculation Type    | Delta_E Target    |
|---------------------|-------------------|
| Standard            | 1 meV/atom        |
| High precision      | 0.1 meV/atom      |
| Force calculations  | 0.01 meV/atom     |
| Phonons             | 0.01 meV/atom     |

---

## 15. Input File Templates

### Minimal SCF Input

```
&CONTROL
  calculation = 'scf'
  prefix = 'silicon'
  outdir = './tmp'
  pseudo_dir = './pseudo'
/
&SYSTEM
  ibrav = 2
  celldm(1) = 10.26
  nat = 2
  ntyp = 1
  ecutwfc = 50.0
/
&ELECTRONS
  conv_thr = 1.0d-8
/
ATOMIC_SPECIES
  Si 28.086 Si.upf
ATOMIC_POSITIONS {crystal}
  Si 0.00 0.00 0.00
  Si 0.25 0.25 0.25
K_POINTS {automatic}
  8 8 8 0 0 0
```

### Band Structure Input (after SCF)

```
&CONTROL
  calculation = 'bands'
  prefix = 'silicon'
  outdir = './tmp'
  pseudo_dir = './pseudo'
/
&SYSTEM
  ibrav = 2
  celldm(1) = 10.26
  nat = 2
  ntyp = 1
  ecutwfc = 50.0
  nbnd = 12
/
&ELECTRONS
  conv_thr = 1.0d-8
/
ATOMIC_SPECIES
  Si 28.086 Si.upf
ATOMIC_POSITIONS {crystal}
  Si 0.00 0.00 0.00
  Si 0.25 0.25 0.25
K_POINTS {crystal_b}
5
  0.000 0.000 0.000 20  ! G
  0.500 0.000 0.500 20  ! X
  0.500 0.250 0.750 20  ! W
  0.375 0.375 0.750 20  ! K
  0.000 0.000 0.000 20  ! G
```

---

## 16. Quick Formulas

### Energy Calculations

```
Formation Energy:
  E_f = E_compound - sum(n_i * E_element_i) / N_atoms

Energy Above Hull:
  E_hull = E_compound - E_convex_hull

Cohesive Energy:
  E_coh = (sum(n_i * E_atom_i) - E_bulk) / N_atoms
```

### Structural Parameters

```
Tolerance Factor (Perovskites):
  t = (r_A + r_O) / [sqrt(2) * (r_B + r_O)]
  Stable perovskite: 0.8 < t < 1.0

Goldschmidt Factor:
  t > 1.0: Hexagonal stacking
  0.9 < t < 1.0: Cubic perovskite
  t < 0.9: Orthorhombic distortion
```

---

## Quick Command Cheat Sheet

```
+------------------+----------------------------------+
| Task             | Command                          |
+------------------+----------------------------------+
| Run SCF          | pw.x -in scf.in > scf.out        |
| Run parallel     | mpirun -np N pw.x -in X > Y      |
| Band structure   | bands.x -in bands_pp.in > X      |
| DOS              | dos.x -in dos.in > dos.out       |
| PDOS             | projwfc.x -in X > Y              |
| Phonons          | ph.x -in ph.in > ph.out          |
| Check output     | grep "!" scf.out                 |
| Check forces     | grep "Total force" relax.out     |
| Check pressure   | grep "P=" vc-relax.out           |
+------------------+----------------------------------+
```

---

*Quick Reference Card for DFT Calculations with Quantum ESPRESSO*
*Version 1.0 - Printable Reference*
