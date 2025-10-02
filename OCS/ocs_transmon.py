"""
Offset-Charge-Sensitive (OCS) Transmon Analysis

This module implements simulation and analysis tools for OCS transmons,
based on the physics described in:
- Serniak et al., PRL (2019): https://arxiv.org/pdf/1903.00113
- Additional context from https://arxiv.org/pdf/2405.17192

The OCS transmon is in the intermediate regime between Cooper-pair box 
(E_J/E_C ≈ 1) and standard transmon (E_J/E_C ≳ 50), where charge 
dispersion is measurable.
"""

import numpy as np
import matplotlib.pyplot as plt
from scipy.constants import h, k, e, eV
from scipy.linalg import eigh


class OCS:
    """
    Offset-Charge-Sensitive Transmon Simulator
    
    This class simulates OCS transmons by diagonalizing the Cooper-pair
    box Hamiltonian and computing dispersive shifts for charge-parity
    readout.
    """
    
    # Physical constants (in SI unless specified)
    PLANCK_EV_S = h / eV  # Planck constant in eV·s
    KB_EV_K = k / eV  # Boltzmann constant in eV/K
    ELECTRON_CHARGE = e  # Elementary charge in Coulombs
    
    # Aluminum material properties
    DOS_AL = 1.72e10  # Density of states [1/(μm³·eV)]
    FERMI_AL = 11.6  # Fermi level [eV]
    TC_AL = 1.2  # Critical temperature [K] (AlMn: 0.1, Al: 1.2)
    DELTA_AL = 1.89e-4  # Superconducting gap [eV]
    
    def __init__(self, e_josephson, e_charging, temperature=0.02,
                 normal_resistance=27e3, delta_left=None, 
                 delta_right=None):
        """
        Initialize OCS transmon simulator
        
        Parameters
        ----------
        e_josephson : float
            Josephson energy [eV]
        e_charging : float
            Charging energy [eV]
        temperature : float, optional
            Operating temperature [K], default 0.02 K
        normal_resistance : float, optional
            Normal state resistance [Ω], default 27 kΩ
        delta_left : float, optional
            Left superconducting gap [eV], defaults to DELTA_AL
        delta_right : float, optional
            Right superconducting gap [eV], defaults to DELTA_AL
        """
        self.e_j = e_josephson
        self.e_c = e_charging
        self.temperature = temperature
        self.r_n = normal_resistance
        self.delta_l = delta_left if delta_left is not None else self.DELTA_AL
        self.delta_r = (delta_right if delta_right is not None 
                       else self.DELTA_AL)
        
        # Computed quantities
        self.ej_ec_ratio = self.e_j / self.e_c
        self.curly_n = (self.DOS_AL * 
                       np.sqrt(2 * np.pi * self.delta_l * 
                               self.KB_EV_K * self.temperature))
        
    @classmethod
    def from_capacitance(cls, total_capacitance, delta, 
                        normal_resistance, **kwargs):
        """
        Create OCS instance from capacitance and resistance
        
        Parameters
        ----------
        total_capacitance : float
            Total capacitance Cᵨ = Cⱼ + Cᵍ + Cₛₕᵤₙₜ [F]
        delta : float
            Superconducting gap [eV]
        normal_resistance : float
            Normal state resistance [Ω]
        **kwargs
            Additional arguments passed to __init__
            
        Returns
        -------
        OCS
            Configured OCS instance
        """
        e_c = cls.ELECTRON_CHARGE / (2 * total_capacitance)  # eV
        e_j = (cls.PLANCK_EV_S * delta / 
              (8 * cls.ELECTRON_CHARGE * normal_resistance))  # eV
        return cls(e_j, e_c, normal_resistance=normal_resistance, 
                  **kwargs)
    
    def build_hamiltonian(self, offset_charge, charge_cutoff=18):
        """
        Construct Cooper-pair box Hamiltonian in charge basis
        
        The Hamiltonian is H = 4Eᴄ(n̂ - nᵍ)² - (Eⱼ/2)(|n⟩⟨n+1| + h.c.)
        where n̂ is the Cooper pair number operator and nᵍ is the 
        offset charge.
        
        Parameters
        ----------
        offset_charge : float
            Dimensionless offset charge nᵍ = CᵍVᵍ/(2e)
        charge_cutoff : int, optional
            Charge basis cutoff, default 18
            
        Returns
        -------
        h : ndarray
            Hamiltonian matrix [eV]
        """
        n_dim = 2 * charge_cutoff + 1
        h = np.zeros((n_dim, n_dim))
        
        # Charge states from -n to +n
        charge_states = np.arange(-charge_cutoff, charge_cutoff + 1)
        
        # Diagonal: charging energy term 4Eᴄ(n - nᵍ)²
        h[np.arange(n_dim), np.arange(n_dim)] = (
            4.0 * self.e_c * (charge_states - offset_charge) ** 2
        )
        
        # Off-diagonal: Josephson tunneling -Eⱼ/2
        h[np.arange(n_dim - 1), np.arange(1, n_dim)] = -self.e_j / 2
        h[np.arange(1, n_dim), np.arange(n_dim - 1)] = -self.e_j / 2
        
        return h
    
    def solve_eigensystem(self, offset_charge, charge_cutoff=18):
        """
        Diagonalize Hamiltonian to find energy eigenstates
        
        Parameters
        ----------
        offset_charge : float
            Dimensionless offset charge nᵍ
        charge_cutoff : int, optional
            Charge basis cutoff, default 18
            
        Returns
        -------
        eigenvalues : ndarray
            Energy eigenvalues [eV], sorted ascending
        eigenvectors : ndarray
            Energy eigenvectors as columns, sorted by eigenvalue
        """
        h = self.build_hamiltonian(offset_charge, charge_cutoff)
        eigenvalues, eigenvectors = eigh(h)
        return eigenvalues, eigenvectors
    
    def solve_system(self, offset_charges, num_levels=4, 
                    charge_cutoff=18):
        """
        Solve for even and odd parity energy levels
        
        Even parity: extra Cooper pair on island (nᵍ)
        Odd parity: extra unpaired electron (nᵍ + 0.5)
        
        Parameters
        ----------
        offset_charges : array_like
            Array of offset charge values to compute
        num_levels : int, optional
            Number of energy levels to return, default 4
        charge_cutoff : int, optional
            Charge basis cutoff, default 18
            
        Returns
        -------
        energies_even : ndarray
            Even parity energies [eV], shape (len(u), num_levels)
        energies_odd : ndarray
            Odd parity energies [eV], shape (len(u), num_levels)
        energy_diff : ndarray
            Parity splitting energy [eV], shape (len(u),)
        """
        offset_charges = np.atleast_1d(offset_charges)
        num_points = len(offset_charges)
        
        energies_even = np.zeros((num_points, num_levels))
        energies_odd = np.zeros((num_points, num_levels))
        
        # Small offset to avoid exact degeneracy
        epsilon = 1e-4
        
        for i, u in enumerate(offset_charges):
            # Even parity: integer charge
            evals_even, _ = self.solve_eigensystem(u + epsilon, 
                                                   charge_cutoff)
            energies_even[i, :] = evals_even[:num_levels]
            
            # Odd parity: half-integer charge
            evals_odd, _ = self.solve_eigensystem(u + 0.5 + epsilon, 
                                                  charge_cutoff)
            energies_odd[i, :] = evals_odd[:num_levels]
        
        # Parity splitting includes superconducting gap difference
        energy_diff = (energies_odd[:, 0] - energies_even[:, 0] + 
                      self.delta_l - self.delta_r)
        
        return energies_even, energies_odd, energy_diff
    
    def compute_dispersive_matrix(self, offset_charge, coupling_g, 
                                  resonator_freq, num_levels=6, 
                                  charge_cutoff=30):
        """
        Compute dispersive shift and matrix elements
        
        The dispersive shift is:
        χᵢ = g² ∑ⱼ≠ᵢ [2ωᵢⱼ |⟨j,p|n̂|i,p⟩|² / (ωᵢⱼ² - ωᵣ²)]
        
        Parameters
        ----------
        offset_charge : float
            Dimensionless offset charge nᵍ
        coupling_g : float
            Transmon-resonator coupling strength [Hz]
        resonator_freq : float
            Resonator frequency [Hz]
        num_levels : int, optional
            Number of levels for matrix elements, default 6
        charge_cutoff : int, optional
            Charge basis cutoff (higher for accuracy), default 30
            
        Returns
        -------
        matrix_elements : ndarray
            Matrix elements |⟨j|n̂|0⟩|², shape (num_levels, num_levels)
        chi : ndarray
            Dispersive shift χᵢ [Hz], shape (num_levels,)
        """
        # Solve eigensystem
        eigenvalues, eigenvectors = self.solve_eigensystem(
            offset_charge, charge_cutoff
        )
        
        # Build number operator n̂ = ∑ₙ (n - nᵍ)|n⟩⟨n|
        n_dim = 2 * charge_cutoff + 1
        charge_states = np.arange(-charge_cutoff, charge_cutoff + 1)
        number_operator = np.diag(charge_states - offset_charge)
        
        # Transform to energy eigenbasis
        total_states = len(eigenvalues)
        matrix_elements = np.zeros((num_levels, num_levels))
        chi = np.zeros(num_levels)
        
        for i in range(num_levels):
            for j in range(total_states):
                if i != j:
                    # Transition frequency
                    omega_ij = (eigenvalues[i] - eigenvalues[j]) / (
                        self.PLANCK_EV_S
                    )  # Hz
                    
                    # Matrix element |⟨j|n̂|i⟩|²
                    mat_elem_sq = np.abs(
                        eigenvectors[:, j].conj() @ number_operator @ 
                        eigenvectors[:, i]
                    ) ** 2
                    
                    # Store matrix elements for first num_levels states
                    if j < num_levels:
                        matrix_elements[i, j] = mat_elem_sq
                    
                    # Dispersive shift contribution
                    chi_contrib = (2.0 * omega_ij * mat_elem_sq / 
                                  (omega_ij ** 2 - resonator_freq ** 2))
                    chi[i] += chi_contrib
        
        # Scale by coupling strength squared
        chi *= coupling_g ** 2
        
        return matrix_elements, chi
    
    def plot_energy_levels(self, offset_charges=None, num_levels=4, 
                          figsize=(10, 7.5)):
        """
        Plot energy level diagram vs offset charge
        
        Shows both even parity (solid) and odd parity (dashed) levels.
        
        Parameters
        ----------
        offset_charges : array_like, optional
            Offset charge values, default linspace(0, 1, 500)
        num_levels : int, optional
            Number of levels to plot, default 4
        figsize : tuple, optional
            Figure size (width, height), default (10, 7.5)
            
        Returns
        -------
        fig : matplotlib.figure.Figure
            Figure handle
        ax : matplotlib.axes.Axes
            Axes handle
        """
        if offset_charges is None:
            offset_charges = np.linspace(0, 1, 500)
        
        energies_even, energies_odd, energy_diff = (
            self.solve_system(offset_charges, num_levels)
        )
        
        # Convert to frequencies (GHz) relative to ground state
        freq_even = ((energies_even - energies_even[:, [0]]) / 
                    self.PLANCK_EV_S / 1e9)
        freq_odd = ((energies_odd - energies_odd[:, [0]]) / 
                   self.PLANCK_EV_S / 1e9)
        
        # Print key frequencies
        print(f"Eⱼ/Eᴄ = {self.ej_ec_ratio:.2f}")
        print(f"E₀₁ = {freq_odd[0, 1]:.3f} GHz")
        print(f"E₀₂ = {freq_odd[0, 2]:.3f} GHz")
        print(f"E₀₃ = {freq_odd[0, 3]:.3f} GHz")
        print(f"ΔE/Eᴄ = {energy_diff[0] / self.e_c:.4f}")
        
        fig, ax = plt.subplots(figsize=figsize)
        
        # Plot even parity (solid lines)
        for j in range(num_levels):
            ax.plot(offset_charges, freq_even[:, j], linewidth=2, 
                   label=f'|{j},e⟩' if j < 2 else None)
        
        # Plot odd parity (dashed lines)
        for j in range(num_levels):
            ax.plot(offset_charges, freq_odd[:, j], '--', linewidth=2,
                   label=f'|{j},o⟩' if j < 2 else None)
        
        ax.set_xlabel(r'Offset Charge [$C_g V_g / 2e$]', fontsize=14)
        ax.set_ylabel(r'$f_{0j}$ [GHz]', fontsize=14)
        ax.set_title(f'$E_J / E_C = {self.ej_ec_ratio:.1f}$', 
                    fontsize=14)
        ax.tick_params(labelsize=12)
        ax.minorticks_on()
        if num_levels <= 4:
            ax.legend(fontsize=12)
        ax.grid(alpha=0.3)
        
        plt.tight_layout()
        return fig, ax
    
    def plot_matrix_elements(self, offset_charges=None, coupling_g=150e6,
                            resonator_freq=7.0e9, num_levels=6, 
                            figsize=(10, 7.5)):
        """
        Plot charge matrix elements vs offset charge
        
        Shows |⟨j,o|n̂|0,o⟩| for transitions from ground state.
        
        Parameters
        ----------
        offset_charges : array_like, optional
            Offset charge values, default linspace(0, 0.5, 250)
        coupling_g : float, optional
            Coupling strength [Hz], default 150 MHz
        resonator_freq : float, optional
            Resonator frequency [Hz], default 7 GHz
        num_levels : int, optional
            Number of levels, default 6
        figsize : tuple, optional
            Figure size, default (10, 7.5)
            
        Returns
        -------
        fig, ax : matplotlib figure and axes
        """
        if offset_charges is None:
            offset_charges = np.linspace(0, 0.5, 250)
        
        num_points = len(offset_charges)
        matrix_elems = np.zeros((num_points, num_levels))
        
        for i, u in enumerate(offset_charges):
            mat, _ = self.compute_dispersive_matrix(
                u + 0.5, coupling_g, resonator_freq, num_levels
            )
            matrix_elems[i, :] = mat[0, :]  # From ground state
        
        fig, ax = plt.subplots(figsize=figsize)
        
        # Plot only transitions (j>0)
        for j in range(1, num_levels):
            ax.semilogy(offset_charges, matrix_elems[:, j], 
                       linewidth=2, label=f'j={j}')
        
        ax.set_ylim([1e-5, 2e0])
        ax.set_xlabel(r'Offset Charge [$C_g V_g / 2e$]', fontsize=14)
        ax.set_ylabel(r'$|\langle j,o|\hat{n}|0,o\rangle|$', 
                     fontsize=14)
        ax.tick_params(labelsize=12)
        ax.minorticks_on()
        ax.legend(fontsize=12)
        ax.grid(alpha=0.3, which='both')
        
        plt.tight_layout()
        return fig, ax
    
    def plot_dispersive_shift(self, offset_charges=None, 
                             coupling_g=150e6, resonator_freq=7.0e9,
                             num_levels=6, figsize=(10, 7.5)):
        """
        Plot dispersive shift χ vs offset charge
        
        Shows charge-parity-dependent shift for ground and first 
        excited states.
        
        Parameters
        ----------
        offset_charges : array_like, optional
            Offset charge values, default linspace(0, 0.5, 250)
        coupling_g : float, optional
            Coupling strength [Hz], default 150 MHz
        resonator_freq : float, optional
            Resonator frequency [Hz], default 7 GHz
        num_levels : int, optional
            Number of levels, default 6
        figsize : tuple, optional
            Figure size, default (10, 7.5)
            
        Returns
        -------
        fig, ax : matplotlib figure and axes
        """
        if offset_charges is None:
            offset_charges = np.linspace(0, 0.5, 250)
        
        num_points = len(offset_charges)
        chi_vals = np.zeros((num_points, num_levels))
        
        for i, u in enumerate(offset_charges):
            _, chi = self.compute_dispersive_matrix(
                u + 0.5, coupling_g, resonator_freq, num_levels
            )
            chi_vals[i, :] = chi
        
        fig, ax = plt.subplots(figsize=figsize)
        
        # Plot first half and mirror for even/odd comparison
        mid = len(offset_charges) // 2
        quarter = len(offset_charges) // 4
        
        # First quarter - odd parity
        for j in range(2):
            ax.plot(offset_charges[:quarter+1], 
                   chi_vals[:quarter+1, j] / 1e6, 
                   linewidth=2, label=f'|{j},o⟩')
        
        # Third quarter mirrored - even parity
        idx_start = mid
        idx_end = mid + quarter + 1
        u_mirrored = offset_charges[idx_start:idx_end] - (
            offset_charges[idx_start]
        )
        for j in range(2):
            ax.plot(u_mirrored, chi_vals[idx_start:idx_end, j] / 1e6,
                   linewidth=2, label=f'|{j},e⟩')
        
        ax.set_xlim([0, 0.25])
        ax.set_ylim([-20, 20])
        ax.set_xlabel(r'Offset Charge [$C_g V_g / 2e$]', fontsize=14)
        ax.set_ylabel(
            r'$\chi_{i,p}$ [MHz]', 
            fontsize=14
        )
        ax.tick_params(labelsize=12)
        ax.minorticks_on()
        ax.legend(fontsize=12)
        ax.grid(alpha=0.3)
        
        plt.tight_layout()
        return fig, ax
    
    def plot_parity_shift_vs_frequency(self, freq_range=None, 
                                      coupling_g=150e6, num_levels=6,
                                      figsize=(10, 7.5)):
        """
        Plot parity-dependent dispersive shift vs resonator frequency
        
        Compares χ₀ at two different offset charges to show maximum
        parity contrast.
        
        Parameters
        ----------
        freq_range : array_like, optional
            Resonator frequencies [Hz], auto-computed if None
        coupling_g : float, optional
            Coupling strength [Hz], default 150 MHz
        num_levels : int, optional
            Number of levels, default 6
        figsize : tuple, optional
            Figure size, default (10, 7.5)
            
        Returns
        -------
        fig, ax : matplotlib figure and axes
        """
        # Auto-determine frequency range if not provided
        if freq_range is None:
            _, energies_odd, _ = self.solve_system([0, 0.5], 4)
            freq_min = 7.0e9  # Start at 7 GHz
            freq_odd_3 = ((energies_odd[0, 3] - energies_odd[0, 0]) / 
                         self.PLANCK_EV_S)  # f_03
            freq_max = freq_odd_3 + 0.2e9  # Add 200 MHz
            freq_range = np.arange(freq_min, freq_max, 1e6)
        
        chi_diff = np.zeros(len(freq_range))
        
        for i, f_r in enumerate(freq_range):
            # Chi at offset charge 0
            _, chi_1 = self.compute_dispersive_matrix(
                0.5, coupling_g, f_r, num_levels
            )
            # Chi at offset charge 0.5
            _, chi_2 = self.compute_dispersive_matrix(
                1.0, coupling_g, f_r, num_levels
            )
            chi_diff[i] = chi_1[0] - chi_2[0]
        
        fig, ax = plt.subplots(figsize=figsize)
        ax.semilogy(freq_range / 1e9, np.abs(chi_diff) / 1e6, 
                   linewidth=2)
        ax.set_xlabel('Resonator Frequency [GHz]', fontsize=14)
        ax.set_ylabel(r'$|\Delta\chi_0|$ [MHz]', fontsize=14)
        ax.tick_params(labelsize=12)
        ax.minorticks_on()
        ax.grid(alpha=0.3, which='both')
        
        # Print summary
        _, energies_odd, _ = self.solve_system([0], 4)
        freq_01 = ((energies_odd[0, 1] - energies_odd[0, 0]) / 
                  self.PLANCK_EV_S)
        freq_02 = ((energies_odd[0, 2] - energies_odd[0, 0]) / 
                  self.PLANCK_EV_S)
        anharmonicity = freq_02 - 2 * freq_01
        
        # Chi at resonator for qubit state readout (simple estimate)
        if hasattr(self, '_resonator_freq'):
            f_r_use = self._resonator_freq
        else:
            f_r_use = 7.0e9
        chi_resonator = (coupling_g ** 2 * anharmonicity / 
                        np.abs(f_r_use - freq_01) / 
                        (np.abs(f_r_use - freq_01) + anharmonicity))
        
        print(f"\nResonator frequency: {f_r_use / 1e9:.2f} GHz")
        print(f"χ (resonator state shift): {chi_resonator / 1e6:.3f} MHz")
        print(f"Δχ (parity shift): {chi_diff[0] / 1e6:.3f} MHz")
        
        plt.tight_layout()
        return fig, ax
    
    def plot_all(self, offset_charges=None, coupling_g=150e6, 
                resonator_freq=7.0e9, num_levels=6):
        """
        Generate all standard plots for OCS transmon analysis
        
        Creates four plots:
        1. Energy level diagram
        2. Matrix elements
        3. Dispersive shift vs offset charge
        4. Parity shift vs resonator frequency
        
        Parameters
        ----------
        offset_charges : array_like, optional
            Offset charge values
        coupling_g : float, optional
            Coupling strength [Hz], default 150 MHz
        resonator_freq : float, optional
            Resonator frequency [Hz], default 7 GHz
        num_levels : int, optional
            Number of levels, default 6
            
        Returns
        -------
        figs : list
            List of figure handles
        """
        self._resonator_freq = resonator_freq  # Store for later use
        
        print("=" * 60)
        print(f"OCS Transmon Parameters:")
        print(f"Eⱼ = {self.e_j / self.KB_EV_K:.3f} K·kᴮ = "
              f"{self.e_j / self.PLANCK_EV_S / 1e9:.3f} GHz")
        print(f"Eᴄ = {self.e_c / self.KB_EV_K:.4f} K·kᴮ = "
              f"{self.e_c / self.PLANCK_EV_S / 1e9:.4f} GHz")
        print(f"Eⱼ/Eᴄ = {self.ej_ec_ratio:.2f}")
        print(f"Rₙ = {self.r_n / 1e3:.1f} kΩ")
        print(f"Resonator frequency = {resonator_freq / 1e9:.2f} GHz")
        print(f"Coupling g = {coupling_g / 1e6:.1f} MHz")
        print(f"FWHM = {resonator_freq / 10000 / 1e6:.2f} MHz")
        print("=" * 60)
        
        figs = []
        
        # Figure 1: Energy levels
        print("\n[1/4] Plotting energy levels...")
        fig1, _ = self.plot_energy_levels(offset_charges, num_levels=4)
        figs.append(fig1)
        
        # Figure 2: Matrix elements
        print("[2/4] Plotting matrix elements...")
        fig2, _ = self.plot_matrix_elements(
            None, coupling_g, resonator_freq, num_levels
        )
        figs.append(fig2)
        
        # Figure 3: Dispersive shift
        print("[3/4] Plotting dispersive shift...")
        fig3, _ = self.plot_dispersive_shift(
            None, coupling_g, resonator_freq, num_levels
        )
        figs.append(fig3)
        
        # Figure 4: Parity shift vs frequency
        print("[4/4] Plotting parity shift vs frequency...")
        fig4, _ = self.plot_parity_shift_vs_frequency(
            None, coupling_g, num_levels
        )
        figs.append(fig4)
        
        print("\nAll plots complete!")
        return figs


def main():
    """Example usage of OCS class"""
    # Example from the MATLAB script (WashU parameters)
    kb_ev_k = OCS.KB_EV_K
    
    ej_ec_ratio = 12
    e_j = 0.4 * kb_ev_k  # eV
    e_c = e_j / ej_ec_ratio  # eV
    
    # Create OCS instance
    ocs = OCS(
        e_josephson=e_j,
        e_charging=e_c,
        temperature=0.02,  # K
        normal_resistance=27e3  # Ω
    )
    
    # Generate all plots
    figs = ocs.plot_all(
        coupling_g=150e6,  # 150 MHz
        resonator_freq=7.0e9,  # 7 GHz
        num_levels=6
    )
    
    plt.show()
    
    return ocs, figs


if __name__ == "__main__":
    ocs, figs = main()

