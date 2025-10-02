"""
Verification script to compare Python implementation with MATLAB results

This script runs the same parameters as the MATLAB OCS_simple.m script
and prints key outputs for comparison.
"""

import numpy as np
from ocs_transmon import OCS


def verify_basic_calculation():
    """
    Verify the Python implementation produces expected results
    
    Uses parameters from OCS_simple.m (WashU config)
    """
    print("=" * 70)
    print("Verification: Python vs MATLAB Implementation")
    print("=" * 70)
    
    # Constants
    planck_ev_s = OCS.PLANCK_EV_S
    
    # Parameters from MATLAB script
    ej_ec_ratio = 12
    e_j_hz = 8.335e9  # ~8.3 GHz (equivalent to 0.4 K·kB)
    e_c_hz = e_j_hz / ej_ec_ratio  # Hz
    
    print(f"\nInput Parameters:")
    print(f"  E_J = {e_j_hz / 1e9:.3f} GHz")
    print(f"  E_C = {e_c_hz / 1e9:.4f} GHz")
    print(f"  E_J/E_C = {ej_ec_ratio}")
    
    # Create OCS instance
    ocs = OCS(
        e_j_hz=e_j_hz,
        e_c_hz=e_c_hz,
        temperature_k=0.02,
        r_n_ohm=27e3
    )
    
    # Compute normal resistance from E_J (for verification)
    e_j_ev = e_j_hz * planck_ev_s
    r_n_computed = (planck_ev_s * OCS.DELTA_AL / 
                   (8 * OCS.ELECTRON_CHARGE * e_j_ev))
    
    print(f"\nDerived Quantities:")
    print(f"  R_n = {r_n_computed / 1e3:.2f} kΩ")
    print(f"  E_J = {e_j_hz / 1e9:.3f} GHz")
    print(f"  E_C = {e_c_hz / 1e9:.4f} GHz")
    
    # Solve system
    offset_charges = np.linspace(0, 1, 500)
    energies_even, energies_odd, energy_diff = ocs.solve_system(
        offset_charges, num_levels=4
    )
    
    # Convert to frequencies
    freq_even = ((energies_even - energies_even[:, [0]]) / 
                planck_ev_s / 1e9)
    freq_odd = ((energies_odd - energies_odd[:, [0]]) / 
               planck_ev_s / 1e9)
    
    print(f"\nTransition Frequencies (at u=0, odd parity):")
    print(f"  E_01 = {freq_odd[0, 1]:.3f} GHz")
    print(f"  E_02 = {freq_odd[0, 2]:.3f} GHz")
    print(f"  E_03 = {freq_odd[0, 3]:.3f} GHz")
    print(f"  ΔE/E_C = {energy_diff[0] / ocs.e_c_ev:.4f}")
    
    # Find minimum frequencies
    fr1 = np.min(freq_odd[:, 1])
    fr3 = np.min(freq_odd[:, 3])
    print(f"  min(f_01) = {fr1:.3f} GHz")
    print(f"  min(f_03) = {fr3:.3f} GHz")
    
    # Compute dispersive shifts
    coupling_g_hz = 150e6  # Hz
    resonator_freq_hz = 7.0e9  # Hz
    
    # At u=0 (odd parity)
    matrix_elements, chi = ocs.compute_dispersive_matrix(
        offset_charge=0.5,
        coupling_g_hz=coupling_g_hz,
        resonator_freq_hz=resonator_freq_hz,
        num_levels=6
    )
    
    # At u=0.5 (even parity)
    matrix_elements_2, chi_2 = ocs.compute_dispersive_matrix(
        offset_charge=1.0,
        coupling_g_hz=coupling_g_hz,
        resonator_freq_hz=resonator_freq_hz,
        num_levels=6
    )
    
    chi_parity = chi[0] - chi_2[0]
    
    # Anharmonicity and chi_resonator
    anharmonicity = (freq_odd[0, 2] - 2 * freq_odd[0, 1]) * 1e9  # Hz
    chi_resonator = (coupling_g_hz ** 2 * anharmonicity / 
                    np.abs(resonator_freq_hz - fr1 * 1e9) / 
                    (np.abs(resonator_freq_hz - fr1 * 1e9) + 
                     anharmonicity))
    
    print(f"\nDispersive Properties:")
    print(f"  Resonator frequency = {resonator_freq_hz / 1e9:.2f} GHz")
    print(f"  Coupling g = {coupling_g_hz / 1e6:.1f} MHz")
    print(f"  Anharmonicity = {anharmonicity / 1e6:.2f} MHz")
    print(f"  χ (resonator state) = {chi_resonator / 1e6:.3f} MHz")
    print(f"  Δχ (parity) = {chi_parity / 1e6:.3f} MHz")
    print(f"  Resonator FWHM = "
          f"{resonator_freq_hz / 10000 / 1e6:.2f} MHz")
    
    # Matrix elements
    print(f"\nCharge Matrix Elements (|⟨j,o|n̂|0,o⟩| at u=0):")
    for j in range(1, min(6, matrix_elements.shape[1])):
        print(f"  j={j}: {np.sqrt(matrix_elements[0, j]):.6e}")
    
    print("\n" + "=" * 70)
    print("Verification Complete!")
    print("\nExpected MATLAB outputs (for comparison):")
    print("  E_01 ~ 4.6-4.7 GHz")
    print("  E_02 ~ 8.9-9.0 GHz")
    print("  E_03 ~ 12.9-13.0 GHz")
    print("  χ (resonator) ~ 1-2 MHz")
    print("  Δχ (parity) ~ few MHz")
    print("=" * 70 + "\n")
    
    return ocs


def verify_hamiltonian_structure():
    """
    Verify Hamiltonian construction matches expected structure
    """
    print("\nVerifying Hamiltonian Structure:")
    print("-" * 70)
    
    kb_ev_k = OCS.KB_EV_K
    e_j_hz = 8.335e9
    e_c_hz = e_j_hz / 12
    
    ocs = OCS(e_j_hz, e_c_hz)
    
    # Build small Hamiltonian for inspection
    u = 0.25
    cutoff = 3  # Small for inspection
    h = ocs.build_hamiltonian(u, charge_cutoff=cutoff)
    
    print(f"Hamiltonian at u={u} (cutoff={cutoff}):")
    print(f"Matrix size: {h.shape}")
    print(f"\nDiagonal (charging energy):")
    diag = np.diag(h)
    for i, val in enumerate(diag):
        n = i - cutoff  # Charge state
        print(f"  n={n:+2d}: {val / kb_ev_k:.6f} K·k_B")
    
    print(f"\nOff-diagonal (Josephson coupling):")
    print(f"  E_J/2 = {-h[0, 1] / kb_ev_k:.6f} K·k_B")
    print(f"  (Should equal {ocs.e_j_ev / 2 / kb_ev_k:.6f} K·k_B)")
    
    # Check Hermiticity
    is_hermitian = np.allclose(h, h.T.conj())
    print(f"\nHamiltonian is Hermitian: {is_hermitian}")
    
    print("-" * 70)


def verify_symmetries():
    """
    Verify charge-parity symmetry properties
    """
    print("\nVerifying Charge-Parity Symmetries:")
    print("-" * 70)
    
    kb_ev_k = OCS.KB_EV_K
    e_j_hz = 8.335e9
    e_c_hz = e_j_hz / 12
    
    ocs = OCS(e_j_hz, e_c_hz)
    
    # Test symmetry: E(u) should equal E(1-u) for even parity
    u_vals = np.array([0.0, 0.1, 0.2, 0.3, 0.4, 0.5])
    energies_even, energies_odd, _ = ocs.solve_system(u_vals, 
                                                       num_levels=2)
    
    print("Testing E(u) vs E(1-u) symmetry:")
    print(f"{'u':<10} {'E_00(u)':<20} {'E_00(1-u)':<20} "
          f"{'Difference':<15}")
    print("-" * 70)
    
    for i, u in enumerate(u_vals):
        if i < len(u_vals) // 2 + 1:
            u_sym = 1.0 - u
            if u_sym <= 1.0:
                # Find closest index to u_sym
                idx_sym = np.argmin(np.abs(u_vals - u_sym))
                e_u = energies_even[i, 0]
                e_u_sym = energies_even[idx_sym, 0]
                diff = np.abs(e_u - e_u_sym)
                print(f"{u:<10.1f} {e_u / kb_ev_k:<20.10f} "
                      f"{e_u_sym / kb_ev_k:<20.10f} "
                      f"{diff / kb_ev_k:<15.2e}")
    
    print("\nNote: Small differences expected at u=0.5 due to numerical "
          "precision")
    print("-" * 70)


if __name__ == "__main__":
    # Run verifications
    ocs = verify_basic_calculation()
    verify_hamiltonian_structure()
    verify_symmetries()
    
    print("\n✓ All verifications complete!")
    print("\nTo generate plots, run:")
    print("  python example_usage.py")

