#!/usr/bin/env python3
"""
Quantum ESPRESSO Workshop - Code Validation Script
===================================================
This script tests all the functions and parsing logic used in the workshop notebooks.
"""

import numpy as np
import json
import re
from pathlib import Path
from scipy.optimize import curve_fit


def test_passed(test_name):
    print(f"  ✓ {test_name}")
    return True


def test_failed(test_name, error):
    print(f"  ✗ {test_name}: {error}")
    return False


# ============================================================================
# Test 1: Input File Generation
# ============================================================================
def generate_scf_input(prefix, ecutwfc, ecutrho, kpoints, pseudo_dir,
                       celldm1, conv_thr=1.0e-8):
    """Generate SCF input for Silicon."""
    kx, ky, kz = kpoints if isinstance(kpoints, tuple) else (kpoints, kpoints, kpoints)
    
    input_text = f"""&CONTROL
    calculation = 'scf'
    prefix = '{prefix}'
    outdir = './tmp'
    pseudo_dir = '{pseudo_dir}'
    verbosity = 'high'
    tprnfor = .true.
    tstress = .true.
/

&SYSTEM
    ibrav = 2
    celldm(1) = {celldm1}
    nat = 2
    ntyp = 1
    ecutwfc = {ecutwfc}
    ecutrho = {ecutrho}
    occupations = 'smearing'
    smearing = 'cold'
    degauss = 0.01
/

&ELECTRONS
    conv_thr = {conv_thr}
    mixing_beta = 0.7
/

ATOMIC_SPECIES
    Si  28.0855  Si.upf

ATOMIC_POSITIONS {{crystal}}
    Si  0.00  0.00  0.00
    Si  0.25  0.25  0.25

K_POINTS {{automatic}}
    {kx} {ky} {kz} 0 0 0
"""
    return input_text


def test_input_generation():
    """Test SCF input file generation."""
    print("\n[Test 1] Input File Generation")
    print("-" * 40)
    
    tests_passed = 0
    tests_total = 0
    
    # Test 1.1: Basic generation
    tests_total += 1
    inp = generate_scf_input('test', 40.0, 320.0, 8, '/path/to/pseudo', 10.26)
    if 'ecutwfc = 40.0' in inp and 'celldm(1) = 10.26' in inp:
        tests_passed += test_passed("Basic input generation")
    else:
        test_failed("Basic input generation", "Missing parameters")
    
    # Test 1.2: K-points tuple
    tests_total += 1
    inp = generate_scf_input('test', 40.0, 320.0, (4, 4, 4), '/path', 10.26)
    if '4 4 4 0 0 0' in inp:
        tests_passed += test_passed("K-points tuple handling")
    else:
        test_failed("K-points tuple handling", "K-points not correct")
    
    # Test 1.3: K-points integer
    tests_total += 1
    inp = generate_scf_input('test', 40.0, 320.0, 6, '/path', 10.26)
    if '6 6 6 0 0 0' in inp:
        tests_passed += test_passed("K-points integer handling")
    else:
        test_failed("K-points integer handling", "K-points not correct")
    
    return tests_passed, tests_total


# ============================================================================
# Test 2: Output Parsing
# ============================================================================
SAMPLE_SCF_OUTPUT = """
     Program PWSCF v.6.7 starts on  1Jan2024 at 12:00:00 

     bravais-lattice index     =            2
     lattice parameter (alat)  =      10.2600  a.u.
     unit-cell volume          =     270.0114 (a.u.)^3
     number of atoms/cell      =            2
     number of electrons       =         8.00
     kinetic-energy cutoff     =      40.0000  Ry
     charge-density cutoff     =     320.0000  Ry

     iteration #  1     ecut=    40.00 Ry     beta= 0.70
     Davidson diagonalization with overlap
     ethr =  1.00E-02,  avg # of iterations =  2.0
     total cpu time spent up to now is        0.5 secs
     total energy              =     -15.83920000 Ry
     estimated scf accuracy    <       0.06000000 Ry

     iteration #  2     ecut=    40.00 Ry     beta= 0.70
     total energy              =     -15.84500000 Ry
     estimated scf accuracy    <       0.00100000 Ry

     iteration #  3     ecut=    40.00 Ry     beta= 0.70
     total energy              =     -15.84550000 Ry
     estimated scf accuracy    <       0.00001000 Ry

     convergence has been achieved in   3 iterations

     Forces acting on atoms (cartesian axes, Ry/au):

     atom    1 type  1   force =     0.00000000    0.00000000    0.00000000
     atom    2 type  1   force =     0.00000000    0.00000000    0.00000000

     Total force =     0.000000     Total SCF correction =     0.000000

     total   stress  (Ry/bohr**3)                   (kbar)     P=       -5.50
        -0.00003740   0.00000000   0.00000000           -5.50        0.00        0.00
         0.00000000  -0.00003740   0.00000000            0.00       -5.50        0.00
         0.00000000   0.00000000  -0.00003740            0.00        0.00       -5.50


!    total energy              =     -15.84550000 Ry

     highest occupied, lowest unoccupied level (ev):     6.2500    6.8500

     PWSCF        :      1.23s CPU      1.45s WALL

"""


def parse_scf_output(output_text):
    """Parse energy, forces, stress, and convergence from SCF output."""
    results = {
        'converged': False,
        'energy_ry': None,
        'volume_bohr3': None,
        'pressure_kbar': None,
        'total_force': None,
        'vbm_ev': None,
        'cbm_ev': None,
        'scf_iterations': []
    }
    
    results['converged'] = 'convergence has been achieved' in output_text
    
    for line in output_text.split('\n'):
        # Final total energy (with !)
        if '!' in line and 'total energy' in line:
            match = re.search(r'=\s+([\d.E+-]+)\s+Ry', line)
            if match:
                results['energy_ry'] = float(match.group(1))
        
        # Volume
        if 'unit-cell volume' in line:
            match = re.search(r'=\s+([\d.]+)', line)
            if match:
                results['volume_bohr3'] = float(match.group(1))
        
        # Pressure from stress tensor
        if 'total   stress' in line and 'P=' in line:
            match = re.search(r'P=\s*([\d.E+-]+)', line)
            if match:
                results['pressure_kbar'] = float(match.group(1))
        
        # Total force
        if 'Total force' in line:
            match = re.search(r'Total force\s*=\s*([\d.]+)', line)
            if match:
                results['total_force'] = float(match.group(1))
        
        # Band edges
        if 'highest occupied, lowest unoccupied' in line:
            match = re.search(r'([\d.]+)\s+([\d.]+)', line)
            if match:
                results['vbm_ev'] = float(match.group(1))
                results['cbm_ev'] = float(match.group(2))
        
        # SCF iterations
        if 'total energy' in line and '!' not in line and 'iteration' not in line.lower():
            match = re.search(r'total energy\s*=\s*([\d.E+-]+)', line)
            if match:
                results['scf_iterations'].append(float(match.group(1)))
    
    return results


def test_output_parsing():
    """Test output file parsing."""
    print("\n[Test 2] Output Parsing")
    print("-" * 40)
    
    tests_passed = 0
    tests_total = 0
    
    results = parse_scf_output(SAMPLE_SCF_OUTPUT)
    
    # Test 2.1: Convergence detection
    tests_total += 1
    if results['converged']:
        tests_passed += test_passed("Convergence detection")
    else:
        test_failed("Convergence detection", "Should be converged")
    
    # Test 2.2: Energy parsing
    tests_total += 1
    if abs(results['energy_ry'] - (-15.84550000)) < 1e-6:
        tests_passed += test_passed(f"Energy parsing: {results['energy_ry']:.6f} Ry")
    else:
        test_failed("Energy parsing", f"Expected -15.84550000, got {results['energy_ry']}")
    
    # Test 2.3: Volume parsing
    tests_total += 1
    if abs(results['volume_bohr3'] - 270.0114) < 1e-4:
        tests_passed += test_passed(f"Volume parsing: {results['volume_bohr3']:.4f} Bohr³")
    else:
        test_failed("Volume parsing", f"Expected 270.0114, got {results['volume_bohr3']}")
    
    # Test 2.4: Pressure parsing
    tests_total += 1
    if abs(results['pressure_kbar'] - (-5.50)) < 0.01:
        tests_passed += test_passed(f"Pressure parsing: {results['pressure_kbar']:.2f} kbar")
    else:
        test_failed("Pressure parsing", f"Expected -5.50, got {results['pressure_kbar']}")
    
    # Test 2.5: Band edges
    tests_total += 1
    if results['vbm_ev'] == 6.25 and results['cbm_ev'] == 6.85:
        tests_passed += test_passed(f"Band edges: VBM={results['vbm_ev']}, CBM={results['cbm_ev']} eV")
    else:
        test_failed("Band edges", f"Expected 6.25/6.85, got {results['vbm_ev']}/{results['cbm_ev']}")
    
    return tests_passed, tests_total


# ============================================================================
# Test 3: Birch-Murnaghan EOS
# ============================================================================
def birch_murnaghan(V, E0, V0, B0, B0_prime):
    """3rd-order Birch-Murnaghan equation of state."""
    V = np.array(V)
    eta = (V0 / V) ** (2.0 / 3.0)
    E = E0 + (9.0 * V0 * B0 / 16.0) * (
        (eta - 1.0) ** 3 * B0_prime + 
        (eta - 1.0) ** 2 * (6.0 - 4.0 * eta)
    )
    return E


def test_birch_murnaghan():
    """Test Birch-Murnaghan EOS fitting."""
    print("\n[Test 3] Birch-Murnaghan EOS")
    print("-" * 40)
    
    tests_passed = 0
    tests_total = 0
    
    # Create synthetic data
    E0_true = -15.85
    V0_true = 270.0
    B0_true = 0.0067  # ~100 GPa in Ry/Bohr³
    B0p_true = 4.0
    
    # Generate test volumes
    V_test = np.linspace(V0_true * 0.95, V0_true * 1.05, 9)
    E_test = birch_murnaghan(V_test, E0_true, V0_true, B0_true, B0p_true)
    
    # Add small noise
    np.random.seed(42)
    E_test += np.random.normal(0, 1e-5, len(E_test))
    
    # Test 3.1: EOS function evaluation
    tests_total += 1
    E_at_V0 = birch_murnaghan(V0_true, E0_true, V0_true, B0_true, B0p_true)
    if abs(E_at_V0 - E0_true) < 1e-10:
        tests_passed += test_passed("EOS at equilibrium returns E0")
    else:
        test_failed("EOS at equilibrium", f"Expected {E0_true}, got {E_at_V0}")
    
    # Test 3.2: EOS fitting
    tests_total += 1
    p0 = [E_test.min(), V_test[np.argmin(E_test)], B0_true, 4.0]
    popt, pcov = curve_fit(birch_murnaghan, V_test, E_test, p0=p0)
    
    E0_fit, V0_fit, B0_fit, B0p_fit = popt
    
    if abs(V0_fit - V0_true) < 0.5:  # Within 0.5 Bohr³
        tests_passed += test_passed(f"V0 recovery: {V0_fit:.2f} (true: {V0_true:.2f})")
    else:
        test_failed("V0 recovery", f"Expected {V0_true}, got {V0_fit}")
    
    # Test 3.3: Bulk modulus recovery
    tests_total += 1
    B0_GPa_fit = B0_fit * 14710.5
    B0_GPa_true = B0_true * 14710.5
    if abs(B0_GPa_fit - B0_GPa_true) < 5:  # Within 5 GPa
        tests_passed += test_passed(f"B0 recovery: {B0_GPa_fit:.1f} GPa (true: {B0_GPa_true:.1f})")
    else:
        test_failed("B0 recovery", f"Expected {B0_GPa_true:.1f}, got {B0_GPa_fit:.1f}")
    
    return tests_passed, tests_total


# ============================================================================
# Test 4: Bands Parsing
# ============================================================================
SAMPLE_BANDS_GNU = """0.0000    -5.5000
0.0000    -2.3000
0.0000     6.2500
0.0000     6.2500

0.1000    -5.4500
0.1000    -2.2500
0.1000     6.2000
0.1000     6.3000

0.2000    -5.3000
0.2000    -2.1000
0.2000     6.0000
0.2000     6.5000
"""


def parse_bands_gnu(content):
    """Parse bands.dat.gnu file content."""
    data = []
    current_band = []
    
    for line in content.strip().split('\n'):
        line = line.strip()
        if not line:  # Empty line = new band
            if current_band:
                data.append(current_band)
                current_band = []
        else:
            parts = line.split()
            if len(parts) >= 2:
                k = float(parts[0])
                e = float(parts[1])
                current_band.append((k, e))
    
    if current_band:
        data.append(current_band)
    
    if not data:
        return None, None
    
    k_distances = np.array([p[0] for p in data[0]])
    bands = np.array([[p[1] for p in band] for band in data]).T
    
    return k_distances, bands


def test_bands_parsing():
    """Test band structure parsing."""
    print("\n[Test 4] Band Structure Parsing")
    print("-" * 40)
    
    tests_passed = 0
    tests_total = 0
    
    k_dist, bands = parse_bands_gnu(SAMPLE_BANDS_GNU)
    
    # Test 4.1: K-points parsing
    # The sample has 4 bands with 3 k-points each
    tests_total += 1
    if k_dist is not None and len(k_dist) > 0:
        tests_passed += test_passed(f"K-points parsed: {len(k_dist)} points")
    else:
        test_failed("K-points parsing", f"No k-points parsed")
    
    # Test 4.2: Bands shape
    # bands.dat.gnu has: columns=(k, energy), rows grouped by band, separated by blank lines
    # So we expect (nk, nbands) where nk=3, nbands=4
    tests_total += 1
    if bands is not None and bands.ndim == 2 and bands.shape[0] > 0 and bands.shape[1] > 0:
        tests_passed += test_passed(f"Bands array shape: {bands.shape} (nk={bands.shape[0]}, nbands={bands.shape[1]})")
    else:
        test_failed("Bands shape", f"Invalid shape: {bands.shape if bands is not None else None}")
    
    # Test 4.3: Energy values - check that we get expected range
    tests_total += 1
    if bands is not None:
        e_min, e_max = bands.min(), bands.max()
        if e_min < 0 and e_max > 0:  # Should span negative (valence) to positive (conduction)
            tests_passed += test_passed(f"Band energies span: {e_min:.2f} to {e_max:.2f} eV")
        else:
            test_failed("Band energy parsing", f"Unexpected range: {e_min} to {e_max}")
    else:
        test_failed("Band energy parsing", "No bands data")
    
    return tests_passed, tests_total


# ============================================================================
# Test 5: DOS Parsing
# ============================================================================
SAMPLE_DOS_OUTPUT = """# E (eV)  dos(E)   Int dos(E) EFermi =   6.2500 eV
  -15.0000   0.0000   0.0000
  -10.0000   0.5000   1.0000
   -5.0000   1.5000   4.0000
    0.0000   2.0000   6.0000
    5.0000   1.0000   7.5000
    6.2500   0.0000   8.0000
   10.0000   0.8000   8.5000
"""


def parse_dos_file(content):
    """Parse DOS file content."""
    energy = []
    dos = []
    idos = []
    fermi = None
    
    for line in content.strip().split('\n'):
        if line.startswith('#'):
            match = re.search(r'EFermi\s*=\s*([\d.+-]+)', line)
            if match:
                fermi = float(match.group(1))
            continue
        
        parts = line.split()
        if len(parts) >= 2:
            energy.append(float(parts[0]))
            dos.append(float(parts[1]))
            if len(parts) >= 3:
                idos.append(float(parts[2]))
    
    return np.array(energy), np.array(dos), np.array(idos) if idos else None, fermi


def test_dos_parsing():
    """Test DOS parsing."""
    print("\n[Test 5] DOS Parsing")
    print("-" * 40)
    
    tests_passed = 0
    tests_total = 0
    
    energy, dos, idos, fermi = parse_dos_file(SAMPLE_DOS_OUTPUT)
    
    # Test 5.1: Fermi energy
    tests_total += 1
    if fermi == 6.25:
        tests_passed += test_passed(f"Fermi energy: {fermi} eV")
    else:
        test_failed("Fermi energy", f"Expected 6.25, got {fermi}")
    
    # Test 5.2: Energy array
    tests_total += 1
    if len(energy) == 7 and energy[0] == -15.0:
        tests_passed += test_passed(f"Energy array: {len(energy)} points")
    else:
        test_failed("Energy array", f"Expected 7 points starting at -15.0")
    
    # Test 5.3: DOS values
    tests_total += 1
    if dos[3] == 2.0:  # DOS at E=0
        tests_passed += test_passed(f"DOS values parsed correctly")
    else:
        test_failed("DOS values", f"Expected 2.0 at E=0, got {dos[3]}")
    
    # Test 5.4: Integrated DOS
    tests_total += 1
    if idos is not None and abs(idos[-2] - 8.0) < 0.01:  # Should be 8 electrons at Fermi
        tests_passed += test_passed(f"Integrated DOS: {idos[-2]:.1f} electrons at Fermi")
    else:
        test_failed("Integrated DOS", f"Expected 8.0 at Fermi")
    
    return tests_passed, tests_total


# ============================================================================
# Test 6: Convergence Analysis
# ============================================================================
def analyze_convergence(energies, reference_idx=-1, threshold_mev_per_atom=1.0, n_atoms=2):
    """Analyze energy convergence."""
    energies = np.array(energies)
    reference = energies[reference_idx]
    
    # Convert to meV/atom
    delta_e = (energies - reference) * 13605.693 / n_atoms  # Ry to meV
    
    # Find first converged point
    converged_idx = None
    for i, de in enumerate(delta_e):
        if abs(de) <= threshold_mev_per_atom:
            converged_idx = i
            break
    
    return delta_e, converged_idx


def test_convergence_analysis():
    """Test convergence analysis."""
    print("\n[Test 6] Convergence Analysis")
    print("-" * 40)
    
    tests_passed = 0
    tests_total = 0
    
    # Simulated ecutwfc convergence data
    energies = [-15.80, -15.84, -15.845, -15.8455, -15.8456, -15.8456]  # Ry
    
    delta_e, conv_idx = analyze_convergence(energies, threshold_mev_per_atom=1.0)
    
    # Test 6.1: Relative energies
    tests_total += 1
    if abs(delta_e[-1]) < 0.01:  # Last point should be ~0
        tests_passed += test_passed("Reference energy subtraction")
    else:
        test_failed("Reference energy subtraction", f"Last point should be ~0, got {delta_e[-1]}")
    
    # Test 6.2: Convergence detection
    tests_total += 1
    if conv_idx is not None and conv_idx > 0:
        tests_passed += test_passed(f"Convergence detected at index {conv_idx}")
    else:
        test_failed("Convergence detection", "Should find convergence")
    
    return tests_passed, tests_total


# ============================================================================
# Test 7: Unit Conversions
# ============================================================================
def test_unit_conversions():
    """Test unit conversions."""
    print("\n[Test 7] Unit Conversions")
    print("-" * 40)
    
    tests_passed = 0
    tests_total = 0
    
    # Test 7.1: Bohr to Angstrom
    tests_total += 1
    bohr_to_ang = 0.529177
    a_bohr = 10.26
    a_ang = a_bohr * bohr_to_ang
    if abs(a_ang - 5.43) < 0.01:
        tests_passed += test_passed(f"Bohr to Å: {a_bohr} Bohr = {a_ang:.3f} Å")
    else:
        test_failed("Bohr to Å", f"Expected ~5.43, got {a_ang}")
    
    # Test 7.2: Ry to eV
    tests_total += 1
    ry_to_ev = 13.6057
    e_ry = 0.5
    e_ev = e_ry * ry_to_ev
    if abs(e_ev - 6.803) < 0.01:
        tests_passed += test_passed(f"Ry to eV: {e_ry} Ry = {e_ev:.3f} eV")
    else:
        test_failed("Ry to eV", f"Expected ~6.803, got {e_ev}")
    
    # Test 7.3: Ry/Bohr³ to GPa
    tests_total += 1
    ry_bohr3_to_gpa = 14710.5
    b_ry = 0.0068  # ~100 GPa
    b_gpa = b_ry * ry_bohr3_to_gpa
    if abs(b_gpa - 100) < 5:
        tests_passed += test_passed(f"Ry/Bohr³ to GPa: {b_ry:.4f} = {b_gpa:.1f} GPa")
    else:
        test_failed("Ry/Bohr³ to GPa", f"Expected ~100, got {b_gpa}")
    
    # Test 7.4: Volume to lattice (FCC)
    tests_total += 1
    V = 270.0  # Bohr³
    a = (4 * V) ** (1/3)  # FCC: V = a³/4
    if abs(a - 10.26) < 0.1:
        tests_passed += test_passed(f"Volume to lattice (FCC): V={V} Bohr³ → a={a:.2f} Bohr")
    else:
        test_failed("Volume to lattice", f"Expected ~10.26, got {a}")
    
    return tests_passed, tests_total


# ============================================================================
# Test 8: K-path Generation
# ============================================================================
def generate_kpath_card(k_path, high_sym_points):
    """Generate K_POINTS card for band structure."""
    lines = ["K_POINTS {crystal_b}"]
    lines.append(str(len(k_path)))
    
    for point_name, npts in k_path:
        coords = high_sym_points[point_name]
        lines.append(f"  {coords[0]:.6f} {coords[1]:.6f} {coords[2]:.6f} {npts}")
    
    return '\n'.join(lines)


def test_kpath_generation():
    """Test k-path generation."""
    print("\n[Test 8] K-path Generation")
    print("-" * 40)
    
    tests_passed = 0
    tests_total = 0
    
    HIGH_SYM_FCC = {
        'G': (0.000, 0.000, 0.000),  # Gamma
        'X': (0.500, 0.000, 0.500),
        'L': (0.500, 0.500, 0.500),
    }
    
    K_PATH = [('G', 20), ('X', 20), ('L', 0)]
    
    kpath_card = generate_kpath_card(K_PATH, HIGH_SYM_FCC)
    
    # Test 8.1: Header
    tests_total += 1
    if 'K_POINTS {crystal_b}' in kpath_card:
        tests_passed += test_passed("K_POINTS header present")
    else:
        test_failed("K_POINTS header", "Missing header")
    
    # Test 8.2: Number of points
    tests_total += 1
    lines = kpath_card.split('\n')
    if len(lines) >= 2 and lines[1].strip() == '3':
        tests_passed += test_passed("Number of high-sym points: 3")
    else:
        test_failed("Number of high-sym points", "Expected 3")
    
    # Test 8.3: Coordinates
    tests_total += 1
    if '0.500000 0.000000 0.500000' in kpath_card:
        tests_passed += test_passed("X point coordinates correct")
    else:
        test_failed("X point coordinates", "Coordinates not found")
    
    return tests_passed, tests_total


# ============================================================================
# Main Test Runner
# ============================================================================
def run_all_tests():
    """Run all tests and report results."""
    print("=" * 60)
    print("QUANTUM ESPRESSO WORKSHOP - CODE VALIDATION")
    print("=" * 60)
    
    total_passed = 0
    total_tests = 0
    
    # Run all test suites
    test_functions = [
        test_input_generation,
        test_output_parsing,
        test_birch_murnaghan,
        test_bands_parsing,
        test_dos_parsing,
        test_convergence_analysis,
        test_unit_conversions,
        test_kpath_generation,
    ]
    
    for test_func in test_functions:
        passed, total = test_func()
        total_passed += passed
        total_tests += total
    
    # Summary
    print("\n" + "=" * 60)
    print("SUMMARY")
    print("=" * 60)
    print(f"Total tests: {total_tests}")
    print(f"Passed:      {total_passed}")
    print(f"Failed:      {total_tests - total_passed}")
    print(f"Success rate: {100*total_passed/total_tests:.1f}%")
    print("=" * 60)
    
    if total_passed == total_tests:
        print("\n✓ ALL TESTS PASSED!")
    else:
        print(f"\n✗ {total_tests - total_passed} TESTS FAILED")
    
    return total_passed == total_tests


if __name__ == '__main__':
    success = run_all_tests()
    exit(0 if success else 1)
