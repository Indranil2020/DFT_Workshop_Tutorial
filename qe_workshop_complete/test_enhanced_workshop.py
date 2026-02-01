#!/usr/bin/env python3
"""
Enhanced Quantum ESPRESSO Workshop - Comprehensive Code Validation
===================================================================
This script tests ALL Python code from the enhanced workshop notebooks.
Tests are organized by notebook and use mock data to avoid actual QE runs.
"""

import numpy as np
import re
import sys
from scipy.optimize import curve_fit


# ============================================================================
# Test Framework
# ============================================================================
class TestResults:
    """Track test results."""
    def __init__(self):
        self.passed = 0
        self.failed = 0
        self.tests = []

    def add_pass(self, name, detail=""):
        self.passed += 1
        self.tests.append(('PASS', name, detail))
        print(f"  [PASS] {name}: {detail}" if detail else f"  [PASS] {name}")

    def add_fail(self, name, detail=""):
        self.failed += 1
        self.tests.append(('FAIL', name, detail))
        print(f"  [FAIL] {name}: {detail}" if detail else f"  [FAIL] {name}")

    @property
    def total(self):
        return self.passed + self.failed


# ============================================================================
# Constants
# ============================================================================
BOHR_TO_ANGSTROM = 0.529177
RY_TO_EV = 13.6057
RY_BOHR3_TO_GPA = 14710.5


# ============================================================================
# Test 1: Unit Conversions (Notebook 00/02)
# ============================================================================
def bohr_to_angstrom(bohr):
    """Convert Bohr to Angstrom."""
    return bohr * BOHR_TO_ANGSTROM


def ry_to_ev(ry):
    """Convert Rydberg to eV."""
    return ry * RY_TO_EV


def ry_bohr3_to_gpa(ry_bohr3):
    """Convert Ry/Bohr^3 to GPa."""
    return ry_bohr3 * RY_BOHR3_TO_GPA


def volume_to_lattice_fcc(V):
    """Convert FCC primitive cell volume to cubic lattice parameter."""
    # FCC primitive cell: V = a^3 / 4
    return (4 * V) ** (1.0 / 3.0)


def test_unit_conversions(results):
    """Test unit conversion functions."""
    print("\n[Test 1] Unit Conversions")
    print("-" * 40)

    # Test 1.1: Bohr to Angstrom
    a_bohr = 10.26
    a_ang = bohr_to_angstrom(a_bohr)
    expected_ang = 5.43
    if abs(a_ang - expected_ang) < 0.01:
        results.add_pass("Bohr to Angstrom", f"{a_bohr} Bohr = {a_ang:.3f} A")
    else:
        results.add_fail("Bohr to Angstrom", f"Expected ~{expected_ang}, got {a_ang}")

    # Test 1.2: Ry to eV
    e_ry = 1.0
    e_ev = ry_to_ev(e_ry)
    if abs(e_ev - 13.606) < 0.01:
        results.add_pass("Ry to eV", f"{e_ry} Ry = {e_ev:.3f} eV")
    else:
        results.add_fail("Ry to eV", f"Expected ~13.606, got {e_ev}")

    # Test 1.3: Ry/Bohr^3 to GPa
    b_ry = 0.0068
    b_gpa = ry_bohr3_to_gpa(b_ry)
    if abs(b_gpa - 100) < 5:
        results.add_pass("Ry/Bohr^3 to GPa", f"{b_ry:.4f} Ry/Bohr^3 = {b_gpa:.1f} GPa")
    else:
        results.add_fail("Ry/Bohr^3 to GPa", f"Expected ~100, got {b_gpa}")

    # Test 1.4: Volume to lattice (FCC)
    V = 270.0  # Bohr^3
    a = volume_to_lattice_fcc(V)
    if abs(a - 10.26) < 0.1:
        results.add_pass("Volume to lattice (FCC)", f"V={V} Bohr^3 -> a={a:.2f} Bohr")
    else:
        results.add_fail("Volume to lattice", f"Expected ~10.26, got {a}")


# ============================================================================
# Test 2: Shannon Ionic Radii (Notebook 02)
# ============================================================================
SHANNON_IONIC_RADII = {
    'Ba': {'+2': {'VI': 1.35, 'XII': 1.61}},
    'Ti': {'+4': {'VI': 0.605}},
    'O': {'-2': {'VI': 1.40}},
    'Sr': {'+2': {'VI': 1.18, 'XII': 1.44}},
    'Pb': {'+2': {'VI': 1.19, 'XII': 1.49}},
    'Zr': {'+4': {'VI': 0.72, 'VIII': 0.84}},
    'Ca': {'+2': {'VI': 1.00, 'XII': 1.34}},
    'La': {'+3': {'VI': 1.032, 'XII': 1.36}},
    'Fe': {'+2': {'VI': 0.78, 'IV': 0.63}, '+3': {'VI': 0.645}},
    'Mn': {'+2': {'VI': 0.83}, '+4': {'VI': 0.53}},
    'Ni': {'+2': {'VI': 0.69}},
}


def get_ionic_radius(element, oxidation_state, coordination='VI'):
    """
    Get Shannon ionic radius for an element.

    Parameters
    ----------
    element : str
        Element symbol (e.g., 'Ba', 'Ti')
    oxidation_state : str
        Oxidation state (e.g., '+2', '-2')
    coordination : str
        Coordination number (e.g., 'VI', 'XII')

    Returns
    -------
    float or None
        Ionic radius in Angstrom
    """
    if element not in SHANNON_IONIC_RADII:
        return None
    if oxidation_state not in SHANNON_IONIC_RADII[element]:
        return None
    if coordination not in SHANNON_IONIC_RADII[element][oxidation_state]:
        return None
    return SHANNON_IONIC_RADII[element][oxidation_state][coordination]


def test_ionic_radii(results):
    """Test Shannon ionic radii lookup."""
    print("\n[Test 2] Shannon Ionic Radii")
    print("-" * 40)

    # Test 2.1: Ba2+ radius
    r_ba = get_ionic_radius('Ba', '+2', 'VI')
    if r_ba is not None and abs(r_ba - 1.35) < 0.01:
        results.add_pass("Ba2+ radius (VI)", f"{r_ba} A")
    else:
        results.add_fail("Ba2+ radius", f"Expected 1.35, got {r_ba}")

    # Test 2.2: Ti4+ radius
    r_ti = get_ionic_radius('Ti', '+4', 'VI')
    if r_ti is not None and abs(r_ti - 0.605) < 0.01:
        results.add_pass("Ti4+ radius (VI)", f"{r_ti} A")
    else:
        results.add_fail("Ti4+ radius", f"Expected 0.605, got {r_ti}")

    # Test 2.3: O2- radius
    r_o = get_ionic_radius('O', '-2', 'VI')
    if r_o is not None and abs(r_o - 1.40) < 0.01:
        results.add_pass("O2- radius (VI)", f"{r_o} A")
    else:
        results.add_fail("O2- radius", f"Expected 1.40, got {r_o}")

    # Test 2.4: XII coordination
    r_ba_xii = get_ionic_radius('Ba', '+2', 'XII')
    if r_ba_xii is not None and abs(r_ba_xii - 1.61) < 0.01:
        results.add_pass("Ba2+ radius (XII)", f"{r_ba_xii} A")
    else:
        results.add_fail("Ba2+ radius (XII)", f"Expected 1.61, got {r_ba_xii}")

    # Test 2.5: Unknown element returns None
    r_unknown = get_ionic_radius('Unobtanium', '+5')
    if r_unknown is None:
        results.add_pass("Unknown element returns None", "Correctly handled")
    else:
        results.add_fail("Unknown element", f"Should return None, got {r_unknown}")


# ============================================================================
# Test 3: Charge Neutrality Check (Notebook 02)
# ============================================================================
def check_charge_neutrality(composition):
    """
    Check if a compound is charge neutral.

    Parameters
    ----------
    composition : list of tuples
        Each tuple: (element, oxidation_state, count)
        Example: [('Ba', +2, 1), ('Ti', +4, 1), ('O', -2, 3)]

    Returns
    -------
    tuple
        (is_neutral, total_charge)
    """
    total_charge = 0
    for element, oxidation, count in composition:
        total_charge += oxidation * count
    return (total_charge == 0, total_charge)


def test_charge_neutrality(results):
    """Test charge neutrality check."""
    print("\n[Test 3] Charge Neutrality Check")
    print("-" * 40)

    # Test 3.1: BaTiO3 (should be neutral)
    bto = [('Ba', +2, 1), ('Ti', +4, 1), ('O', -2, 3)]
    is_neutral, charge = check_charge_neutrality(bto)
    if is_neutral and charge == 0:
        results.add_pass("BaTiO3 neutrality", "Ba2+ + Ti4+ + 3*O2- = 0")
    else:
        results.add_fail("BaTiO3 neutrality", f"Total charge = {charge}")

    # Test 3.2: SrTiO3 (should be neutral)
    sto = [('Sr', +2, 1), ('Ti', +4, 1), ('O', -2, 3)]
    is_neutral, charge = check_charge_neutrality(sto)
    if is_neutral:
        results.add_pass("SrTiO3 neutrality", "Sr2+ + Ti4+ + 3*O2- = 0")
    else:
        results.add_fail("SrTiO3 neutrality", f"Total charge = {charge}")

    # Test 3.3: Invalid compound (should not be neutral)
    invalid = [('Ba', +2, 1), ('Ti', +4, 1), ('O', -2, 2)]  # BaTiO2 - wrong
    is_neutral, charge = check_charge_neutrality(invalid)
    if not is_neutral and charge == +2:
        results.add_pass("Invalid compound detection", f"BaTiO2 has charge +{charge}")
    else:
        results.add_fail("Invalid compound", f"Should be non-neutral, charge={charge}")

    # Test 3.4: La2O3 (should be neutral)
    la2o3 = [('La', +3, 2), ('O', -2, 3)]
    is_neutral, charge = check_charge_neutrality(la2o3)
    if is_neutral:
        results.add_pass("La2O3 neutrality", "2*La3+ + 3*O2- = 0")
    else:
        results.add_fail("La2O3 neutrality", f"Total charge = {charge}")


# ============================================================================
# Test 4: Lattice Parameter Estimation (Notebook 02)
# ============================================================================
def estimate_lattice_ionic_radii(r_A, r_B, r_O=1.40):
    """
    Estimate perovskite lattice parameter from ionic radii.
    For ideal cubic perovskite ABO3: a = r_A + r_O (A-O distance * sqrt(2))
    More accurate: a = sqrt(2) * (r_A + r_O)
    """
    return np.sqrt(2) * (r_A + r_O)


def estimate_lattice_isostructural(a_ref, V_ref, V_new):
    """
    Estimate lattice parameter by isostructural scaling.
    a_new = a_ref * (V_new/V_ref)^(1/3)
    """
    return a_ref * (V_new / V_ref) ** (1.0 / 3.0)


def test_lattice_estimation(results):
    """Test lattice parameter estimation methods."""
    print("\n[Test 4] Lattice Parameter Estimation")
    print("-" * 40)

    # Test 4.1: Ionic radii sum for BaTiO3
    r_Ba = 1.61  # XII coordination
    r_O = 1.40
    a_estimated = estimate_lattice_ionic_radii(r_Ba, 0.605, r_O)
    # BaTiO3 exp. a ~ 4.00 A
    if 3.5 < a_estimated < 4.5:
        results.add_pass("Ionic radii sum method", f"BaTiO3 a ~ {a_estimated:.2f} A")
    else:
        results.add_fail("Ionic radii sum", f"Got {a_estimated}, expected ~4.0")

    # Test 4.2: Isostructural scaling
    a_ref = 3.905  # SrTiO3 lattice parameter
    V_ref = a_ref ** 3
    V_new = 4.0 ** 3  # Slightly larger volume
    a_new = estimate_lattice_isostructural(a_ref, V_ref, V_new)
    expected = 4.0
    if abs(a_new - expected) < 0.01:
        results.add_pass("Isostructural scaling", f"a_new = {a_new:.3f} A")
    else:
        results.add_fail("Isostructural scaling", f"Expected {expected}, got {a_new}")


# ============================================================================
# Test 5: Input File Generation (Notebooks 03-06)
# ============================================================================
def generate_scf_input(prefix, ecutwfc, ecutrho, kpoints, pseudo_dir,
                       celldm1=10.26, conv_thr=1.0e-8):
    """Generate SCF input file for Silicon."""
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


def generate_relax_input(prefix, ecutwfc, ecutrho, kpoints, pseudo_dir, celldm1=10.26):
    """Generate relax input file."""
    kx, ky, kz = kpoints if isinstance(kpoints, tuple) else (kpoints, kpoints, kpoints)

    input_text = f"""&CONTROL
    calculation = 'relax'
    prefix = '{prefix}'
    outdir = './tmp'
    pseudo_dir = '{pseudo_dir}'
    forc_conv_thr = 1.0e-4
/

&SYSTEM
    ibrav = 2
    celldm(1) = {celldm1}
    nat = 2
    ntyp = 1
    ecutwfc = {ecutwfc}
    ecutrho = {ecutrho}
/

&ELECTRONS
    conv_thr = 1.0e-8
/

&IONS
    ion_dynamics = 'bfgs'
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


def generate_bands_input(prefix, ecutwfc, ecutrho, pseudo_dir, celldm1, kpath_card, nbnd=12):
    """Generate bands input file."""
    input_text = f"""&CONTROL
    calculation = 'bands'
    prefix = '{prefix}'
    outdir = './tmp'
    pseudo_dir = '{pseudo_dir}'
/

&SYSTEM
    ibrav = 2
    celldm(1) = {celldm1}
    nat = 2
    ntyp = 1
    ecutwfc = {ecutwfc}
    ecutrho = {ecutrho}
    nbnd = {nbnd}
/

&ELECTRONS
    conv_thr = 1.0e-8
/

ATOMIC_SPECIES
    Si  28.0855  Si.upf

ATOMIC_POSITIONS {{crystal}}
    Si  0.00  0.00  0.00
    Si  0.25  0.25  0.25

{kpath_card}
"""
    return input_text


def generate_magnetic_input(prefix, ecutwfc, ecutrho, kpoints, pseudo_dir,
                           celldm1, starting_mag):
    """Generate magnetic SCF input file."""
    kx, ky, kz = kpoints if isinstance(kpoints, tuple) else (kpoints, kpoints, kpoints)

    mag_lines = "\n".join([f"    starting_magnetization({i+1}) = {m}"
                          for i, m in enumerate(starting_mag)])

    input_text = f"""&CONTROL
    calculation = 'scf'
    prefix = '{prefix}'
    outdir = './tmp'
    pseudo_dir = '{pseudo_dir}'
/

&SYSTEM
    ibrav = 2
    celldm(1) = {celldm1}
    nat = 2
    ntyp = 1
    ecutwfc = {ecutwfc}
    ecutrho = {ecutrho}
    nspin = 2
{mag_lines}
/

&ELECTRONS
    conv_thr = 1.0e-8
/

ATOMIC_SPECIES
    Fe  55.845  Fe.upf

ATOMIC_POSITIONS {{crystal}}
    Fe  0.00  0.00  0.00
    Fe  0.50  0.50  0.50

K_POINTS {{automatic}}
    {kx} {ky} {kz} 0 0 0
"""
    return input_text


def test_input_generation(results):
    """Test input file generation."""
    print("\n[Test 5] Input File Generation")
    print("-" * 40)

    # Test 5.1: SCF input generation
    inp = generate_scf_input('test', 40.0, 320.0, (8, 8, 8), '/path/to/pseudo', 10.26)
    has_ecutwfc = 'ecutwfc = 40.0' in inp
    has_celldm = 'celldm(1) = 10.26' in inp
    has_kpts = '8 8 8 0 0 0' in inp
    if has_ecutwfc and has_celldm and has_kpts:
        results.add_pass("SCF input generation", "All parameters present")
    else:
        results.add_fail("SCF input generation", "Missing parameters")

    # Test 5.2: Relax input generation
    inp = generate_relax_input('test', 40.0, 320.0, 6, '/path', 10.26)
    has_relax = "calculation = 'relax'" in inp
    has_ions = '&IONS' in inp
    has_bfgs = 'bfgs' in inp
    if has_relax and has_ions and has_bfgs:
        results.add_pass("Relax input generation", "Contains IONS namelist")
    else:
        results.add_fail("Relax input generation", "Missing IONS namelist")

    # Test 5.3: Bands input generation
    kpath = "K_POINTS {crystal_b}\n3\n0.0 0.0 0.0 20\n0.5 0.0 0.5 20\n0.5 0.5 0.5 0"
    inp = generate_bands_input('test', 40.0, 320.0, '/path', 10.26, kpath, 12)
    has_bands = "calculation = 'bands'" in inp
    has_nbnd = 'nbnd = 12' in inp
    has_kpath = 'crystal_b' in inp
    if has_bands and has_nbnd and has_kpath:
        results.add_pass("Bands input generation", "Contains nbnd and k-path")
    else:
        results.add_fail("Bands input generation", "Missing bands parameters")

    # Test 5.4: Magnetic input generation
    inp = generate_magnetic_input('fe', 60.0, 480.0, 8, '/path', 5.42, [0.5, 0.5])
    has_nspin = 'nspin = 2' in inp
    has_mag = 'starting_magnetization' in inp
    if has_nspin and has_mag:
        results.add_pass("Magnetic input generation", "Contains nspin=2 and starting_mag")
    else:
        results.add_fail("Magnetic input generation", "Missing magnetic parameters")


# ============================================================================
# Test 6: Output Parsing (All Notebooks)
# ============================================================================
MOCK_SCF_OUTPUT = """
     Program PWSCF v.7.2 starts on  1Jan2024 at 12:00:00

     bravais-lattice index     =            2
     lattice parameter (alat)  =      10.2600  a.u.
     unit-cell volume          =     270.0114 (a.u.)^3
     number of atoms/cell      =            2
     number of electrons       =         8.00
     kinetic-energy cutoff     =      40.0000  Ry
     charge-density cutoff     =     320.0000  Ry
     number of k points=     60

     iteration #  1     ecut=    40.00 Ry     beta= 0.70
     Davidson diagonalization with overlap
     ethr =  1.00E-02,  avg # of iterations =  2.0
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

     atom    1 type  1   force =     0.00001234    0.00000000    0.00000000
     atom    2 type  1   force =    -0.00001234    0.00000000    0.00000000

     Total force =     0.000017     Total SCF correction =     0.000000

     total   stress  (Ry/bohr**3)                   (kbar)     P=       -5.50
        -0.00003740   0.00000000   0.00000000           -5.50        0.00        0.00
         0.00000000  -0.00003740   0.00000000            0.00       -5.50        0.00
         0.00000000   0.00000000  -0.00003740            0.00        0.00       -5.50


!    total energy              =     -15.84550000 Ry

     highest occupied, lowest unoccupied level (ev):     6.2500    6.8500

     PWSCF        :      1.23s CPU      1.45s WALL

"""

MOCK_BANDS_GNU = """0.0000    -5.5000
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

0.3000    -5.1000
0.3000    -1.9000
0.3000     5.8000
0.3000     6.7000
"""

MOCK_DOS_OUTPUT = """# E (eV)  dos(E)   Int dos(E) EFermi =   6.2500 eV
  -15.0000   0.0000   0.0000
  -10.0000   0.5000   1.0000
   -5.0000   1.5000   4.0000
    0.0000   2.0000   6.0000
    5.0000   1.0000   7.5000
    6.2500   0.0000   8.0000
   10.0000   0.8000   8.5000
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
        'scf_iterations': [],
        'n_kpoints': None,
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

        # Number of k-points
        if 'number of k points' in line:
            match = re.search(r'=\s*(\d+)', line)
            if match:
                results['n_kpoints'] = int(match.group(1))

    return results


def parse_bands_gnu(content):
    """Parse bands.dat.gnu file content."""
    data = []
    current_band = []

    for line in content.strip().split('\n'):
        line = line.strip()
        if not line:
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


def test_output_parsing(results):
    """Test output file parsing."""
    print("\n[Test 6] Output Parsing")
    print("-" * 40)

    # Test 6.1: SCF convergence detection
    scf = parse_scf_output(MOCK_SCF_OUTPUT)
    if scf['converged']:
        results.add_pass("SCF convergence detection", "Detected convergence")
    else:
        results.add_fail("SCF convergence detection", "Should be converged")

    # Test 6.2: SCF energy parsing
    if scf['energy_ry'] is not None and abs(scf['energy_ry'] - (-15.84550000)) < 1e-6:
        results.add_pass("SCF energy parsing", f"{scf['energy_ry']:.8f} Ry")
    else:
        results.add_fail("SCF energy parsing", f"Expected -15.84550000, got {scf['energy_ry']}")

    # Test 6.3: Volume parsing
    if scf['volume_bohr3'] is not None and abs(scf['volume_bohr3'] - 270.0114) < 0.001:
        results.add_pass("Volume parsing", f"{scf['volume_bohr3']:.4f} Bohr^3")
    else:
        results.add_fail("Volume parsing", f"Expected 270.0114, got {scf['volume_bohr3']}")

    # Test 6.4: Pressure parsing
    if scf['pressure_kbar'] is not None and abs(scf['pressure_kbar'] - (-5.50)) < 0.01:
        results.add_pass("Pressure parsing", f"{scf['pressure_kbar']:.2f} kbar")
    else:
        results.add_fail("Pressure parsing", f"Expected -5.50, got {scf['pressure_kbar']}")

    # Test 6.5: Force parsing
    if scf['total_force'] is not None and abs(scf['total_force'] - 0.000017) < 0.0001:
        results.add_pass("Force parsing", f"{scf['total_force']:.6f} Ry/Bohr")
    else:
        results.add_fail("Force parsing", f"Expected ~0.000017, got {scf['total_force']}")

    # Test 6.6: Band edges parsing
    if scf['vbm_ev'] == 6.25 and scf['cbm_ev'] == 6.85:
        results.add_pass("Band edges parsing", f"VBM={scf['vbm_ev']}, CBM={scf['cbm_ev']} eV")
    else:
        results.add_fail("Band edges parsing", f"VBM={scf['vbm_ev']}, CBM={scf['cbm_ev']}")

    # Test 6.7: Bands parsing
    k_dist, bands = parse_bands_gnu(MOCK_BANDS_GNU)
    if k_dist is not None and len(k_dist) == 4 and bands.shape[1] == 4:
        results.add_pass("Bands parsing", f"{len(k_dist)} k-points, {bands.shape[1]} bands")
    else:
        results.add_fail("Bands parsing", "Incorrect shape")

    # Test 6.8: DOS parsing
    energy, dos, idos, fermi = parse_dos_file(MOCK_DOS_OUTPUT)
    if fermi is not None and abs(fermi - 6.25) < 0.001:
        results.add_pass("DOS Fermi energy", f"{fermi} eV")
    else:
        results.add_fail("DOS Fermi energy", f"Expected 6.25, got {fermi}")

    # Test 6.9: Integrated DOS
    if idos is not None and abs(idos[-2] - 8.0) < 0.01:
        results.add_pass("Integrated DOS at Fermi", f"{idos[-2]:.1f} electrons")
    else:
        results.add_fail("Integrated DOS", f"Expected 8.0 at Fermi")


# ============================================================================
# Test 7: Convergence Analysis (Notebook 04)
# ============================================================================
def analyze_convergence(energies, reference_idx=-1, threshold_mev_per_atom=1.0, n_atoms=2):
    """
    Analyze energy convergence.

    Parameters
    ----------
    energies : array
        Total energies in Ry
    reference_idx : int
        Index of reference energy (default: last)
    threshold_mev_per_atom : float
        Convergence threshold in meV/atom
    n_atoms : int
        Number of atoms

    Returns
    -------
    delta_e : array
        Energy differences in meV/atom
    converged_idx : int or None
        Index of first converged point
    """
    energies = np.array(energies)
    reference = energies[reference_idx]

    # Convert to meV/atom
    delta_e = (energies - reference) * RY_TO_EV * 1000 / n_atoms

    # Find first converged point
    converged_idx = None
    for i, de in enumerate(delta_e):
        if abs(de) <= threshold_mev_per_atom:
            converged_idx = i
            break

    return delta_e, converged_idx


def test_convergence_analysis(results):
    """Test convergence analysis."""
    print("\n[Test 7] Convergence Analysis")
    print("-" * 40)

    # Test 7.1: Energy convergence with mock data
    # Simulated ecutwfc convergence
    energies = [-15.80, -15.84, -15.845, -15.8455, -15.8456, -15.8456]
    delta_e, conv_idx = analyze_convergence(energies, threshold_mev_per_atom=1.0)

    if abs(delta_e[-1]) < 0.01:
        results.add_pass("Reference energy subtraction", "Last point is ~0 meV/atom")
    else:
        results.add_fail("Reference energy", f"Last point should be ~0, got {delta_e[-1]}")

    # Test 7.2: Convergence detection
    if conv_idx is not None and conv_idx > 0:
        results.add_pass("Energy convergence detection", f"Converged at index {conv_idx}")
    else:
        results.add_fail("Convergence detection", "Should find convergence")

    # Test 7.3: K-point convergence
    kpt_energies = [-15.84, -15.8456, -15.8457, -15.8457]  # Converges quickly
    delta_e_k, conv_idx_k = analyze_convergence(kpt_energies, threshold_mev_per_atom=0.5)
    if conv_idx_k is not None:
        results.add_pass("K-point convergence detection", f"Converged at index {conv_idx_k}")
    else:
        results.add_fail("K-point convergence", "Should converge")


# ============================================================================
# Test 8: Birch-Murnaghan EOS (Notebook 05)
# ============================================================================
def birch_murnaghan(V, E0, V0, B0, B0_prime):
    """
    Third-order Birch-Murnaghan equation of state.

    Parameters
    ----------
    V : array
        Volume
    E0 : float
        Equilibrium energy (Ry)
    V0 : float
        Equilibrium volume (Bohr^3)
    B0 : float
        Bulk modulus (Ry/Bohr^3)
    B0_prime : float
        Pressure derivative of bulk modulus

    Returns
    -------
    E : array
        Energy at given volumes (Ry)
    """
    V = np.array(V)
    eta = (V0 / V) ** (2.0 / 3.0)
    E = E0 + (9.0 * V0 * B0 / 16.0) * (
        (eta - 1.0) ** 3 * B0_prime +
        (eta - 1.0) ** 2 * (6.0 - 4.0 * eta)
    )
    return E


def test_birch_murnaghan(results):
    """Test Birch-Murnaghan EOS."""
    print("\n[Test 8] Birch-Murnaghan EOS")
    print("-" * 40)

    # Known parameters
    E0_true = -15.85
    V0_true = 270.0
    B0_true = 0.0067  # ~100 GPa in Ry/Bohr^3
    B0p_true = 4.0

    # Test 8.1: EOS at equilibrium
    E_at_V0 = birch_murnaghan(V0_true, E0_true, V0_true, B0_true, B0p_true)
    if abs(E_at_V0 - E0_true) < 1e-10:
        results.add_pass("EOS at equilibrium", f"E(V0) = E0 = {E0_true} Ry")
    else:
        results.add_fail("EOS at equilibrium", f"Expected {E0_true}, got {E_at_V0}")

    # Test 8.2: Generate and fit synthetic data
    V_test = np.linspace(V0_true * 0.95, V0_true * 1.05, 9)
    E_test = birch_murnaghan(V_test, E0_true, V0_true, B0_true, B0p_true)

    # Add small noise
    np.random.seed(42)
    E_test += np.random.normal(0, 1e-5, len(E_test))

    # Fit
    p0 = [E_test.min(), V_test[np.argmin(E_test)], B0_true, 4.0]
    popt, pcov = curve_fit(birch_murnaghan, V_test, E_test, p0=p0)
    E0_fit, V0_fit, B0_fit, B0p_fit = popt

    # Test 8.3: V0 recovery
    if abs(V0_fit - V0_true) < 0.5:
        results.add_pass("V0 recovery", f"V0 = {V0_fit:.2f} Bohr^3 (true: {V0_true})")
    else:
        results.add_fail("V0 recovery", f"Expected {V0_true}, got {V0_fit}")

    # Test 8.4: Bulk modulus recovery
    B0_GPa_fit = B0_fit * RY_BOHR3_TO_GPA
    B0_GPa_true = B0_true * RY_BOHR3_TO_GPA
    if abs(B0_GPa_fit - B0_GPa_true) < 5:
        results.add_pass("B0 recovery", f"B0 = {B0_GPa_fit:.1f} GPa (true: {B0_GPa_true:.1f})")
    else:
        results.add_fail("B0 recovery", f"Expected {B0_GPa_true:.1f}, got {B0_GPa_fit:.1f}")

    # Test 8.5: Lattice parameter from V0
    a0_fit = volume_to_lattice_fcc(V0_fit)
    a0_true = volume_to_lattice_fcc(V0_true)
    if abs(a0_fit - a0_true) < 0.02:
        results.add_pass("Lattice from EOS", f"a0 = {a0_fit:.3f} Bohr")
    else:
        results.add_fail("Lattice from EOS", f"Expected {a0_true:.3f}, got {a0_fit:.3f}")


# ============================================================================
# Test 9: Elastic Constants (Notebook 07)
# ============================================================================
def check_born_stability_cubic(C11, C12, C44):
    """
    Check Born stability criteria for cubic crystal.

    Criteria:
    1. C11 - C12 > 0
    2. C11 + 2*C12 > 0
    3. C44 > 0
    """
    crit1 = C11 - C12 > 0
    crit2 = C11 + 2 * C12 > 0
    crit3 = C44 > 0
    return crit1 and crit2 and crit3


def check_born_stability_hexagonal(C11, C12, C13, C33, C44):
    """
    Check Born stability criteria for hexagonal crystal.

    Criteria:
    1. C11 > |C12|
    2. (C11 + C12) * C33 > 2 * C13^2
    3. C44 > 0
    """
    crit1 = C11 > abs(C12)
    crit2 = (C11 + C12) * C33 > 2 * C13 ** 2
    crit3 = C44 > 0
    return crit1 and crit2 and crit3


def bulk_modulus_from_cij_cubic(C11, C12):
    """Calculate bulk modulus from elastic constants (Voigt average for cubic)."""
    return (C11 + 2 * C12) / 3


def bulk_modulus_from_cij_hexagonal(C11, C12, C13, C33):
    """Calculate bulk modulus from elastic constants (Voigt average for hexagonal)."""
    return (2 * C11 + 2 * C12 + 4 * C13 + C33) / 9


def test_elastic_constants(results):
    """Test elastic constant analysis."""
    print("\n[Test 9] Elastic Constants")
    print("-" * 40)

    # Test 9.1: Born stability for silicon (cubic)
    # Typical Si values in GPa
    C11_si = 166.0
    C12_si = 64.0
    C44_si = 80.0
    is_stable = check_born_stability_cubic(C11_si, C12_si, C44_si)
    if is_stable:
        results.add_pass("Born stability (cubic Si)", "All criteria satisfied")
    else:
        results.add_fail("Born stability (cubic)", "Should be stable")

    # Test 9.2: Unstable cubic crystal
    C11_bad = 50.0
    C12_bad = 60.0  # C11 - C12 < 0, unstable!
    C44_bad = 30.0
    is_stable_bad = check_born_stability_cubic(C11_bad, C12_bad, C44_bad)
    if not is_stable_bad:
        results.add_pass("Unstable cubic detection", "C11-C12<0 detected")
    else:
        results.add_fail("Unstable detection", "Should be unstable")

    # Test 9.3: Born stability for hexagonal
    # Typical Ti values in GPa
    C11_ti = 162.0
    C12_ti = 92.0
    C13_ti = 69.0
    C33_ti = 181.0
    C44_ti = 47.0
    is_stable_hex = check_born_stability_hexagonal(C11_ti, C12_ti, C13_ti, C33_ti, C44_ti)
    if is_stable_hex:
        results.add_pass("Born stability (hexagonal Ti)", "All criteria satisfied")
    else:
        results.add_fail("Born stability (hexagonal)", "Should be stable")

    # Test 9.4: Bulk modulus from Cij (cubic)
    B_si = bulk_modulus_from_cij_cubic(C11_si, C12_si)
    expected_B = 98.0  # Approx. for Si
    if abs(B_si - expected_B) < 5:
        results.add_pass("Bulk modulus from Cij (cubic)", f"B = {B_si:.1f} GPa")
    else:
        results.add_fail("Bulk modulus (cubic)", f"Expected ~{expected_B}, got {B_si}")


# ============================================================================
# Test 10: K-path Generation (Notebook 08)
# ============================================================================
HIGH_SYM_FCC = {
    'G': (0.000, 0.000, 0.000),
    'X': (0.500, 0.000, 0.500),
    'W': (0.500, 0.250, 0.750),
    'K': (0.375, 0.375, 0.750),
    'L': (0.500, 0.500, 0.500),
}

HIGH_SYM_BCC = {
    'G': (0.000, 0.000, 0.000),
    'H': (0.500, -0.500, 0.500),
    'N': (0.000, 0.000, 0.500),
    'P': (0.250, 0.250, 0.250),
}


def generate_kpath_card(k_path, high_sym_points):
    """Generate K_POINTS {crystal_b} card for band structure."""
    lines = ["K_POINTS {crystal_b}"]
    lines.append(str(len(k_path)))

    for point_name, npts in k_path:
        coords = high_sym_points[point_name]
        lines.append(f"  {coords[0]:.6f} {coords[1]:.6f} {coords[2]:.6f} {npts}")

    return '\n'.join(lines)


def test_kpath_generation(results):
    """Test k-path generation."""
    print("\n[Test 10] K-path Generation")
    print("-" * 40)

    # Test 10.1: FCC k-path
    K_PATH_FCC = [('G', 20), ('X', 20), ('L', 0)]
    kpath_card = generate_kpath_card(K_PATH_FCC, HIGH_SYM_FCC)

    has_header = 'K_POINTS {crystal_b}' in kpath_card
    has_npoints = '3' in kpath_card.split('\n')[1]
    has_gamma = '0.000000 0.000000 0.000000' in kpath_card
    has_x = '0.500000 0.000000 0.500000' in kpath_card

    if has_header and has_npoints and has_gamma and has_x:
        results.add_pass("FCC k-path generation", "G-X-L path correct")
    else:
        results.add_fail("FCC k-path", "Missing elements")

    # Test 10.2: BCC k-path
    K_PATH_BCC = [('G', 20), ('H', 20), ('N', 0)]
    kpath_card_bcc = generate_kpath_card(K_PATH_BCC, HIGH_SYM_BCC)

    has_h = '0.500000 -0.500000 0.500000' in kpath_card_bcc
    if has_h:
        results.add_pass("BCC k-path generation", "H point at (0.5,-0.5,0.5)")
    else:
        results.add_fail("BCC k-path", "H point incorrect")

    # Test 10.3: Correct number of points
    lines = kpath_card.split('\n')
    n_sym_points = int(lines[1])
    if n_sym_points == len(K_PATH_FCC):
        results.add_pass("K-path point count", f"{n_sym_points} high-symmetry points")
    else:
        results.add_fail("K-path point count", f"Expected {len(K_PATH_FCC)}, got {n_sym_points}")


# ============================================================================
# Test 11: Effective Mass Calculation (Notebook 08)
# ============================================================================
def fit_parabolic_band(k, E):
    """
    Fit parabolic band E = E0 + hbar^2 k^2 / (2 m*)

    Parameters
    ----------
    k : array
        k-values (1/Bohr)
    E : array
        Energies (eV)

    Returns
    -------
    E0 : float
        Band extremum (eV)
    m_eff : float
        Effective mass in units of electron mass
    """
    # E = E0 + A * k^2
    # A = hbar^2 / (2 * m* * m_e)
    # In atomic units with E in eV: A = 7.62 / m* (approximately)

    def parabola(k, E0, A):
        return E0 + A * k ** 2

    popt, _ = curve_fit(parabola, k, E)
    E0, A = popt

    # Convert A to effective mass
    # hbar^2 / (2 m_e) = 3.81 eV*A^2 = 7.62 eV*Bohr^2
    hbar2_over_2me = 7.62  # eV * Bohr^2
    m_eff = hbar2_over_2me / A if A != 0 else float('inf')

    return E0, m_eff


def test_effective_mass(results):
    """Test effective mass calculation."""
    print("\n[Test 11] Effective Mass Calculation")
    print("-" * 40)

    # Test 11.1: Known effective mass recovery
    # Create parabolic band with known m* = 0.5
    m_star_true = 0.5
    hbar2_over_2me = 7.62
    A_true = hbar2_over_2me / m_star_true
    E0_true = 0.0

    k_test = np.linspace(-0.1, 0.1, 21)
    E_test = E0_true + A_true * k_test ** 2

    E0_fit, m_star_fit = fit_parabolic_band(k_test, E_test)

    if abs(E0_fit - E0_true) < 0.001:
        results.add_pass("Band extremum recovery", f"E0 = {E0_fit:.4f} eV")
    else:
        results.add_fail("Band extremum", f"Expected {E0_true}, got {E0_fit}")

    if abs(m_star_fit - m_star_true) < 0.01:
        results.add_pass("Effective mass recovery", f"m* = {m_star_fit:.3f} m_e")
    else:
        results.add_fail("Effective mass", f"Expected {m_star_true}, got {m_star_fit}")

    # Test 11.2: Electron-like band (positive curvature)
    k_e = np.linspace(-0.05, 0.05, 11)
    E_e = 0.5 + 15.0 * k_e ** 2  # Electron-like
    E0_e, m_e = fit_parabolic_band(k_e, E_e)
    if m_e > 0:
        results.add_pass("Electron-like band", f"m* = {m_e:.3f} m_e > 0")
    else:
        results.add_fail("Electron-like band", "Should have positive m*")


# ============================================================================
# Test 12: Phonon Thermodynamics (Notebook 09)
# ============================================================================
def phonon_heat_capacity(frequencies_thz, T, n_atoms=2):
    """
    Calculate heat capacity from phonon frequencies using quantum harmonic oscillator.

    Parameters
    ----------
    frequencies_thz : array
        Phonon frequencies in THz
    T : float
        Temperature in K
    n_atoms : int
        Number of atoms

    Returns
    -------
    Cv : float
        Heat capacity in J/(mol*K)
    """
    if T <= 0:
        return 0.0

    # Constants
    kB = 1.380649e-23  # J/K
    hbar = 1.054572e-34  # J*s
    NA = 6.02214076e23  # Avogadro

    # Convert THz to rad/s
    omega = np.array(frequencies_thz) * 2 * np.pi * 1e12

    # Remove acoustic modes at gamma (omega ~ 0)
    omega = omega[omega > 1e10]

    if len(omega) == 0:
        return 0.0

    Cv = 0.0
    for w in omega:
        x = hbar * w / (kB * T)
        if x < 50:  # Avoid overflow
            exp_x = np.exp(x)
            Cv += kB * x ** 2 * exp_x / (exp_x - 1) ** 2

    return Cv * NA


def debye_temperature(v_sound, V_per_atom):
    """
    Calculate Debye temperature from sound velocity.

    Parameters
    ----------
    v_sound : float
        Average sound velocity in m/s
    V_per_atom : float
        Volume per atom in m^3

    Returns
    -------
    theta_D : float
        Debye temperature in K
    """
    hbar = 1.054572e-34
    kB = 1.380649e-23

    # omega_D = v_sound * (6 * pi^2 / V)^(1/3)
    omega_D = v_sound * (6 * np.pi ** 2 / V_per_atom) ** (1.0 / 3.0)
    theta_D = hbar * omega_D / kB

    return theta_D


def test_phonon_thermodynamics(results):
    """Test phonon thermodynamics calculations."""
    print("\n[Test 12] Phonon Thermodynamics")
    print("-" * 40)

    # Test 12.1: Heat capacity at high T (classical limit)
    # At high T, Cv -> 3 N kB per atom = 3R per mole
    # For 2 atoms: Cv -> 6R ~ 50 J/(mol*K)
    freqs = [5.0, 10.0, 15.0, 8.0, 12.0, 14.0]  # 6 modes for 2 atoms
    Cv_high = phonon_heat_capacity(freqs, T=3000, n_atoms=2)
    R = 8.314  # J/(mol*K)
    classical_limit = 6 * R  # 6 modes * R

    # At high T, should approach classical limit
    if Cv_high > 0 and Cv_high < classical_limit * 1.2:
        results.add_pass("Heat capacity high-T", f"Cv(3000K) = {Cv_high:.1f} J/(mol*K)")
    else:
        results.add_fail("Heat capacity high-T", f"Got {Cv_high}, expected < {classical_limit}")

    # Test 12.2: Heat capacity at low T (should be small)
    Cv_low = phonon_heat_capacity(freqs, T=10, n_atoms=2)
    if Cv_low >= 0 and Cv_low < Cv_high:
        results.add_pass("Heat capacity low-T", f"Cv(10K) = {Cv_low:.3f} J/(mol*K)")
    else:
        results.add_fail("Heat capacity low-T", "Should be smaller than high-T value")

    # Test 12.3: Debye temperature
    v_sound = 5000  # m/s (typical for Si)
    V_per_atom = 20e-30  # m^3 (typical for Si)
    theta_D = debye_temperature(v_sound, V_per_atom)
    # Si Debye temperature ~ 640 K
    if 300 < theta_D < 1000:
        results.add_pass("Debye temperature", f"theta_D = {theta_D:.0f} K")
    else:
        results.add_fail("Debye temperature", f"Got {theta_D}, expected 300-1000 K")


# ============================================================================
# Main Test Runner
# ============================================================================
def run_all_tests():
    """Run all tests and report results."""
    print("=" * 40)
    print("ENHANCED QE WORKSHOP - CODE VALIDATION")
    print("=" * 40)

    results = TestResults()

    # Run all test suites
    test_unit_conversions(results)
    test_ionic_radii(results)
    test_charge_neutrality(results)
    test_lattice_estimation(results)
    test_input_generation(results)
    test_output_parsing(results)
    test_convergence_analysis(results)
    test_birch_murnaghan(results)
    test_elastic_constants(results)
    test_kpath_generation(results)
    test_effective_mass(results)
    test_phonon_thermodynamics(results)

    # Summary
    print("\n" + "=" * 40)
    print("SUMMARY")
    print("=" * 40)
    print(f"Total: {results.total}")
    print(f"Passed: {results.passed}")
    print(f"Failed: {results.failed}")

    success_rate = 100 * results.passed / results.total if results.total > 0 else 0
    print(f"Success: {success_rate:.1f}%")

    print("=" * 40)

    if results.failed == 0:
        print("\n[PASS] ALL TESTS PASSED!")
        return 0
    else:
        print(f"\n[FAIL] {results.failed} TESTS FAILED")
        return 1


if __name__ == '__main__':
    exit_code = run_all_tests()
    sys.exit(exit_code)
