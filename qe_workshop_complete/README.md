# Quantum ESPRESSO Workshop

A comprehensive workshop for learning DFT calculations with Quantum ESPRESSO.

## Workshop Structure

### Notebooks (in order)

1. **01_Introduction_and_Setup.ipynb** - Environment setup, pseudopotentials, basic concepts
2. **02_SCF_Calculation_Basics.ipynb** - Input/output structure, parameter explanations
3. **03_Ecutwfc_Convergence.ipynb** - Wavefunction cutoff convergence testing
4. **04_Kpoint_Convergence.ipynb** - Brillouin zone sampling convergence
5. **05_Lattice_Optimization.ipynb** - Equation of state, equilibrium lattice
6. **06_Band_Structure.ipynb** - Electronic band structure calculations
7. **07_DOS_Calculation.ipynb** - Density of states and PDOS
8. **08_Summary_and_Exercises.ipynb** - Summary and practice exercises

### Supporting Files

- `QUICK_REFERENCE.md` - Quick reference card with formulas and parameters
- `test_all_code.py` - Validation script for all workshop code
- `pseudopotentials/` - Si pseudopotential file

## Prerequisites

- Python 3.8+
- numpy, scipy, matplotlib
- Quantum ESPRESSO 6.x or 7.x
- Jupyter Notebook/Lab or Google Colab

## Installation

```bash
# Install Python packages
pip install numpy scipy matplotlib jupyter

# Set up QE (Ubuntu/Debian)
sudo apt install quantum-espresso
```

## Usage

1. Start with Notebook 01 to verify your setup
2. Follow notebooks in sequence
3. Complete convergence testing (03-05) before any property calculations
4. Use the Quick Reference card for parameters and formulas

## Key Principle

**Always converge these parameters BEFORE calculating properties:**
1. ecutwfc (wavefunction cutoff)
2. k-points (Brillouin zone sampling)
3. Lattice parameter (equilibrium structure)

## Validation

Run the test script to verify all code functions correctly:
```bash
python test_all_code.py
```

Expected output: "ALL TESTS PASSED!"

## License

This workshop is provided for educational purposes.

## Author

Generated for DFT workshop training.
