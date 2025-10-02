"""
Offset-Charge-Sensitive (OCS) Transmon Analysis

This module implements simulation and analysis tools for OCS transmons,
based on the physics described in:
- Serniak et al., PRA (2019): https://arxiv.org/pdf/1903.00113
- Additional context from https://arxiv.org/pdf/2405.17192

The OCS transmon is in the intermediate regime between Cooper-pair box 
(E_J/E_C ≈ 1) and standard transmon (E_J/E_C ≳ 50), where charge 
dispersion is measurable.
"""

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm
from scipy.constants import h, k, e, eV
from scipy.linalg import eigh
import yaml
from pathlib import Path


class OCS:
    """
    Offset-Charge-Sensitive Transmon Simulator
    
    This class simulates OCS transmons by diagonalizing the Cooper-pair
    box Hamiltonian and computing dispersive shifts for charge-parity
    readout.
    
    Material properties are loaded from materials.yaml, supporting
    multiple superconductors (Al, Hf, Nb, TiN, etc.)
    """
    
    # Physical constants (in SI unless specified)
    PLANCK_EV_S = h / eV  # Planck constant in eV·s
    KB_EV_K = k / eV  # Boltzmann constant in eV/K
    ELECTRON_CHARGE = e  # Elementary charge in Coulombs
    
    # Materials database cache
    _materials_db = None
    _materials_path = Path(__file__).parent / "materials.yaml"
    _style_path = Path(__file__).parent / "ocs.mplstyle"
    
    @classmethod
    def load_materials_database(cls):
        """
        Load materials database from YAML file
        
        Returns
        -------
        dict
            Dictionary containing material properties
        """
        if cls._materials_db is None:
            with open(cls._materials_path, 'r') as f:
                cls._materials_db = yaml.safe_load(f)
        return cls._materials_db
    
    @classmethod
    def list_materials(cls):
        """
        List available materials in the database
        
        Returns
        -------
        list
            List of available material names
        """
        db = cls.load_materials_database()
        return list(db['materials'].keys())
    
    @classmethod
    def get_material_properties(cls, material_name):
        """
        Get properties for a specific material
        
        Parameters
        ----------
        material_name : str
            Name of the material (e.g., 'aluminum', 'hafnium')
            
        Returns
        -------
        dict
            Material properties dictionary
        """
        db = cls.load_materials_database()
        if material_name not in db['materials']:
            available = ', '.join(cls.list_materials())
            raise ValueError(
                f"Material '{material_name}' not found. "
                f"Available: {available}"
            )
        return db['materials'][material_name]['properties']
    
    def __init__(self, 
                 e_j_hz, e_c_hz, temperature_k=0.02, r_n_ohm=27e3, 
                 delta_l_hz=None, delta_r_hz=None,
                 material='aluminum', **material_overrides):
        """
        Initialize OCS transmon simulator
        
        Parameters
        ----------
        e_j_hz : float
            Josephson energy [Hz]
        e_c_hz : float
            Charging energy [Hz]
        temperature_k : float, optional
            Operating temperature [K], default 0.02 K
        r_n_ohm : float, optional
            Normal state resistance [Ω], default 27 kΩ
        delta_l_hz : float, optional
            Left superconducting gap [Hz], defaults to material value
        delta_r_hz : float, optional
            Right superconducting gap [Hz], defaults to material value
        material : str, optional
            Material name from database, default 'aluminum'
            Options: 'aluminum', 'hafnium', 'niobium', etc.
        **material_overrides : dict
            Override material properties (dos, fermi_level, tc, delta)
            Values should be in SI/eV units as specified in YAML
        """
        # Store input parameters in their native units
        self.e_j_hz = e_j_hz
        self.e_c_hz = e_c_hz
        self.temperature_k = temperature_k
        self.r_n_ohm = r_n_ohm
        
        # Load material properties
        self.material_name = material
        mat_props = self.get_material_properties(material).copy()
        
        # Apply any manual overrides
        mat_props.update(material_overrides)
        
        # Store material properties as instance variables
        # Convert to float in case YAML parsed as string
        self.dos = float(mat_props['dos'])  # [1/(μm³·eV)]
        self.fermi_level = float(mat_props['fermi_level'])  # [eV]
        self.tc = float(mat_props['tc'])  # [K]
        self.delta_material = float(mat_props['delta'])  # [eV]
        
        # Convert to eV for internal calculations
        self.e_j_ev = e_j_hz * self.PLANCK_EV_S
        self.e_c_ev = e_c_hz * self.PLANCK_EV_S
        
        # Delta in eV (convert from Hz if provided, else use material)
        if delta_l_hz is not None:
            self.delta_l_ev = delta_l_hz * self.PLANCK_EV_S
        else:
            self.delta_l_ev = self.delta_material
            
        if delta_r_hz is not None:
            self.delta_r_ev = delta_r_hz * self.PLANCK_EV_S
        else:
            self.delta_r_ev = self.delta_material
        
        # Computed quantities
        self.ej_ec_ratio = self.e_j_hz / self.e_c_hz
        self.curly_n = (self.dos * 
                       np.sqrt(2 * np.pi * self.delta_l_ev * 
                               self.KB_EV_K * self.temperature_k))
        
    @classmethod
    def from_capacitance(cls, c_total_f, delta_hz, r_n_ohm, **kwargs):
        """
        Create OCS instance from capacitance and resistance
        
        Parameters
        ----------
        c_total_f : float
            Total capacitance Cᵨ = Cⱼ + Cᵍ + Cₛₕᵤₙₜ [F]
        delta_hz : float
            Superconducting gap [Hz]
        r_n_ohm : float
            Normal state resistance [Ω]
        **kwargs
            Additional arguments passed to __init__
            
        Returns
        -------
        OCS
            Configured OCS instance
        """
        # Compute energies in eV first
        e_c_ev = cls.ELECTRON_CHARGE / (2 * c_total_f)  # eV
        delta_ev = delta_hz * cls.PLANCK_EV_S  # eV
        e_j_ev = (cls.PLANCK_EV_S * delta_ev / 
                 (8 * cls.ELECTRON_CHARGE * r_n_ohm))  # eV
        
        # Convert to Hz for the constructor
        e_c_hz = e_c_ev / cls.PLANCK_EV_S
        e_j_hz = e_j_ev / cls.PLANCK_EV_S
        
        return cls(e_j_hz, e_c_hz, r_n_ohm=r_n_ohm, **kwargs)
    
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
            4.0 * self.e_c_ev * (charge_states - offset_charge) ** 2
        )
        
        # Off-diagonal: Josephson tunneling -Eⱼ/2
        h[np.arange(n_dim - 1), np.arange(1, n_dim)] = -self.e_j_ev / 2
        h[np.arange(1, n_dim), np.arange(n_dim - 1)] = -self.e_j_ev / 2
        
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
                      self.delta_l_ev - self.delta_r_ev)
        
        return energies_even, energies_odd, energy_diff
    
    def compute_dispersive_matrix(self, offset_charge, coupling_g_hz, 
                                  resonator_freq_hz, num_levels=6, 
                                  charge_cutoff=30):
        """
        Compute dispersive shift and matrix elements
        
        The dispersive shift is:
        χᵢ = g² ∑ⱼ≠ᵢ [2ωᵢⱼ |⟨j,p|n̂|i,p⟩|² / (ωᵢⱼ² - ωᵣ²)]
        
        Parameters
        ----------
        offset_charge : float
            Dimensionless offset charge nᵍ
        coupling_g_hz : float
            Transmon-resonator coupling strength [Hz]
        resonator_freq_hz : float
            Resonator frequency [Hz]
        num_levels : int, optional
            Number of levels for matrix elements, default 6
        charge_cutoff : int, optional
            Charge basis cutoff (higher for accuracy), default 30
            
        Returns
        -------
        matrix_elements : ndarray
            Matrix elements |⟨j|n̂|0⟩|², shape (num_levels, num_levels)
        chi_ip : ndarray
            Dispersive shift χᵢ,ₚ [Hz], shape (num_levels,)
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
        chi_ip = np.zeros(num_levels)
        
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
                                  (omega_ij ** 2 - resonator_freq_hz ** 2))
                    chi_ip[i] += chi_contrib
        
        # Scale by coupling strength squared
        chi_ip *= coupling_g_hz ** 2
        
        return matrix_elements, chi_ip
    
    def plot_energy_levels(self, offset_charges=None, num_levels=5, 
                          freq_resonator_hz=None, coupling_g_hz=None,
                          figsize=(4, 3)):
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
            Figure size (width, height), default (4, 3)
            
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
        print(f"ΔE/Eᴄ = {energy_diff[0] / self.e_c_ev:.4f}")
        
        with plt.style.context(self._style_path):
            fig, ax = plt.subplots(figsize=figsize)
            
            # Get colormaps
            cmap_even = cm.get_cmap('Reds')
            cmap_odd = cm.get_cmap('Blues')
            
            # Plot even parity (solid lines) - Reds colormap
            for j in range(num_levels):
                color = cmap_even(0.2 + 0.7 * j / max(num_levels - 1, 1))
                ax.plot(offset_charges, freq_even[:, j], linewidth=2, 
                       color=color,
                       label=f'|{j},e⟩')
            
            # Plot odd parity (dashed lines) - Blues colormap
            for j in range(num_levels):
                color = cmap_odd(0.3 + 0.7 * j / max(num_levels - 1, 1))
                ax.plot(offset_charges, freq_odd[:, j], 
                       linewidth=2, color=color,
                       label=f'|{j},o⟩')
            
            if freq_resonator_hz is not None:
                ax.axhline(freq_resonator_hz / 1e9, 
                          color='black', linestyle=':')

            ax.set_xlim([0, 1])
            ax.set_xlabel(r'Offset Charge [$C_g V_g / 2e$]')
            ax.set_ylabel(r'$f_{0j}$ [GHz]')
            
            # Construct comprehensive title
            title_parts = [
                f'$\\xi={self.ej_ec_ratio:.1f}$',
                f'$E_J={self.e_j_hz/1e9:.2f}$ GHz',
                f'$E_C={self.e_c_hz/1e9:.3f}$ GHz'
            ]
            if coupling_g_hz is not None:
                title_parts.append(f'$g={coupling_g_hz/1e6:.0f}$ MHz')
            title_parts.append(f'$T={self.temperature_k*1e3:.0f}$ mK')
            ax.set_title(', '.join(title_parts), fontsize=7)
            
            ax.minorticks_on()
            if num_levels <= 4:
                ax.legend(loc="best")
            ax.grid(alpha=0.3)
        
        return fig, ax
    
    def plot_matrix_elements(self, offset_charges=None, 
                            coupling_g_hz=150e6, 
                            resonator_freq_hz=7.0e9, num_levels=6, 
                            figsize=(4, 3)):
        """
        Plot charge matrix elements vs offset charge
        
        Shows |⟨j,o|n̂|0,o⟩| for transitions from ground state.
        
        Parameters
        ----------
        offset_charges : array_like, optional
            Offset charge values, default linspace(0, 1, 500)
        coupling_g_hz : float, optional
            Coupling strength [Hz], default 150 MHz
        resonator_freq_hz : float, optional
            Resonator frequency [Hz], default 7 GHz
        num_levels : int, optional
            Number of levels, default 6
        figsize : tuple, optional
            Figure size, default (4, 3)
            
        Returns
        -------
        fig, ax : matplotlib figure and axes
        """
        if offset_charges is None:
            offset_charges = np.linspace(0, 1, 500)
        
        num_points = len(offset_charges)
        matrix_elems = np.zeros((num_points, num_levels))
        
        for i, u in enumerate(offset_charges):
            mat, _ = self.compute_dispersive_matrix(
                u + 0.5, coupling_g_hz, resonator_freq_hz, num_levels
            )
            matrix_elems[i, :] = mat[0, :]  # From ground state
        
        with plt.style.context(self._style_path):
            fig, ax = plt.subplots(figsize=figsize)
            
            # Get colormap for odd parity
            cmap_odd = cm.get_cmap('Blues')
            
            # Plot only transitions (j>0)
            for j in range(1, num_levels):
                color = cmap_odd(0.3 + 0.7 * j / max(num_levels - 1, 1))
                ax.semilogy(offset_charges, matrix_elems[:, j], 
                           linewidth=2, color=color, label=f'j={j}')
            
            ax.set_ylim([1e-5, 2e0])
            ax.set_xlim([0, 1])
            ax.set_xlabel(r'Offset Charge [$C_g V_g / 2e$]')
            ax.set_ylabel(r'$|\langle j,o|\hat{n}|0,o\rangle|^2$')
            
            # Construct comprehensive title
            title_parts = [
                f'$\\xi={self.ej_ec_ratio:.1f}$',
                f'$E_J={self.e_j_hz/1e9:.2f}$ GHz',
                f'$E_C={self.e_c_hz/1e9:.3f}$ GHz',
                f'$g={coupling_g_hz/1e6:.0f}$ MHz',
                f'$T={self.temperature_k*1e3:.0f}$ mK'
            ]
            ax.set_title(', '.join(title_parts), fontsize=7)
            
            ax.minorticks_on()
            ax.legend(loc="best")
            ax.grid(alpha=0.3, which='both')
        
        return fig, ax
    
    def plot_dispersive_shift(self, offset_charges=None, 
                             coupling_g_hz=150e6, 
                             resonator_freq_hz=7.0e9,
                             num_levels=6, figsize=(4, 3)):
        """
        Plot dispersive shift χ vs offset charge
        
        Shows charge-parity-dependent shift for ground and first 
        excited states.
        
        Parameters
        ----------
        offset_charges : array_like, optional
            Offset charge values, default linspace(0, 1, 500)
        coupling_g_hz : float, optional
            Coupling strength [Hz], default 150 MHz
        resonator_freq_hz : float, optional
            Resonator frequency [Hz], default 7 GHz
        num_levels : int, optional
            Number of levels, default 6
        figsize : tuple, optional
            Figure size, default (4, 3)
            
        Returns
        -------
        fig, ax : matplotlib figure and axes
        """
        if offset_charges is None:
            offset_charges = np.linspace(0, 1, 500)
        
        num_points = len(offset_charges)
        chi_vals = np.zeros((num_points, num_levels))
        
        for i, u in enumerate(offset_charges):
            _, chi_ip = self.compute_dispersive_matrix(
                u + 0.5, coupling_g_hz, resonator_freq_hz, num_levels
            )
            chi_vals[i, :] = chi_ip
        
        with plt.style.context(self._style_path):
            fig, ax = plt.subplots(figsize=figsize)
            
            # Get colormap for odd parity
            cmap_odd = cm.get_cmap('Blues')
            
            # Plot full range showing even and odd parity behavior
            for j in range(num_levels-2):
                color = cmap_odd(0.3 + 0.7 * j / max(num_levels - 3, 1))
                ax.plot(offset_charges, chi_vals[:, j] / 1e6, 
                       linewidth=2, color=color, label=f'|{j}⟩')
            
            ax.set_xlim([0, 1])
            ax.set_ylim([-10, 10])
            ax.set_xlabel(r'Offset Charge [$C_g V_g / 2e$]')
            ax.set_ylabel(r'$\chi_{i,p}$ [MHz]')
            
            # Construct comprehensive title
            title_parts = [
                f'$\\xi={self.ej_ec_ratio:.1f}$',
                f'$E_J={self.e_j_hz/1e9:.2f}$ GHz',
                f'$E_C={self.e_c_hz/1e9:.3f}$ GHz',
                f'$g={coupling_g_hz/1e6:.0f}$ MHz',
                f'$T={self.temperature_k*1e3:.0f}$ mK'
            ]
            ax.set_title(', '.join(title_parts), fontsize=7)
            
            ax.minorticks_on()
            ax.legend(loc="best")
            ax.grid(alpha=0.3)
        
        return fig, ax
    
    def plot_parity_shift_vs_frequency(self, freq_range_hz=None, 
                                      coupling_g_hz=150e6, 
                                      num_levels=6,
                                      figsize=(4, 3)):
        """
        Plot parity-dependent dispersive shift vs resonator frequency
        
        Compares χ₀ at two different offset charges to show maximum
        parity contrast.
        
        Parameters
        ----------
        freq_range_hz : array_like, optional
            Resonator frequencies [Hz], auto-computed if None
        coupling_g_hz : float, optional
            Coupling strength [Hz], default 150 MHz
        num_levels : int, optional
            Number of levels, default 6
        figsize : tuple, optional
            Figure size, default (4, 3)
            
        Returns
        -------
        fig, ax : matplotlib figure and axes
        """
        # Auto-determine frequency range if not provided
        if freq_range_hz is None:
            _, energies_odd, _ = self.solve_system([0, 0.5], 4)
            freq_min = 7.0e9  # Start at 7 GHz
            freq_odd_3 = ((energies_odd[0, 3] - energies_odd[0, 0]) / 
                         self.PLANCK_EV_S)  # f_03
            freq_max = freq_odd_3 + 0.2e9  # Add 200 MHz
            freq_range_hz = np.arange(freq_min, freq_max, 1e6)
        
        chi_diff = np.zeros(len(freq_range_hz))
        
        for i, f_r in enumerate(freq_range_hz):
            # Chi at offset charge 0
            _, chi_ip_1 = self.compute_dispersive_matrix(
                0.5, coupling_g_hz, f_r, num_levels
            )
            # Chi at offset charge 0.5
            _, chi_ip_2 = self.compute_dispersive_matrix(
                1.0, coupling_g_hz, f_r, num_levels
            )
            chi_diff[i] = chi_ip_1[0] - chi_ip_2[0]
        
        with plt.style.context(self._style_path):
            fig, ax = plt.subplots(figsize=figsize)
            
            ax.semilogy(freq_range_hz / 1e9, np.abs(chi_diff) / 1e6, 
                       linewidth=2)
            ax.set_xlabel('Resonator Frequency [GHz]')
            ax.set_ylabel(r'$|\Delta\chi_0|$ [MHz]')
            
            # Construct comprehensive title
            title_parts = [
                f'$\\xi={self.ej_ec_ratio:.1f}$',
                f'$E_J={self.e_j_hz/1e9:.2f}$ GHz',
                f'$E_C={self.e_c_hz/1e9:.3f}$ GHz',
                f'$g={coupling_g_hz/1e6:.0f}$ MHz',
                f'$T={self.temperature_k*1e3:.0f}$ mK'
            ]
            ax.set_title(', '.join(title_parts), fontsize=7)
            
            ax.minorticks_on()
            ax.legend(loc="best")
            ax.grid(alpha=0.3, which='both')
        
        # Print summary
        _, energies_odd, _ = self.solve_system([0], 4)
        freq_01 = ((energies_odd[0, 1] - energies_odd[0, 0]) / 
                  self.PLANCK_EV_S)
        freq_02 = ((energies_odd[0, 2] - energies_odd[0, 0]) / 
                  self.PLANCK_EV_S)
        anharmonicity = freq_02 - 2 * freq_01
        
        # Chi at resonator for qubit state readout (simple estimate)
        if hasattr(self, '_resonator_freq_hz'):
            f_r_use = self._resonator_freq_hz
        else:
            f_r_use = 7.0e9
        chi_resonator = (coupling_g_hz ** 2 * anharmonicity / 
                        np.abs(f_r_use - freq_01) / 
                        (np.abs(f_r_use - freq_01) + anharmonicity))
        
        print(f"\nResonator frequency: {f_r_use / 1e9:.2f} GHz")
        print(f"χ (resonator state shift): {chi_resonator / 1e6:.3f} MHz")
        print(f"Δχ (parity shift): {chi_diff[0] / 1e6:.3f} MHz")
        
        return fig, ax
    
    def plot_all(self, offset_charges=None, coupling_g_hz=150e6, 
                resonator_freq_hz=7.0e9, num_levels=5):
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
        coupling_g_hz : float, optional
            Coupling strength [Hz], default 150 MHz
        resonator_freq_hz : float, optional
            Resonator frequency [Hz], default 7 GHz
        num_levels : int, optional
            Number of levels, default 6
            
        Returns
        -------
        figs : list
            List of figure handles
        """
        self._resonator_freq_hz = resonator_freq_hz  # Store for later use
        
        print("=" * 60)
        print(f"OCS Transmon Parameters:")
        print(f"Material: {self.material_name} (Tc = {self.tc:.3f} K)")
        print(f"Eⱼ = {self.e_j_ev / self.KB_EV_K:.3f} K·kᴮ = "
              f"{self.e_j_hz / 1e9:.3f} GHz")
        print(f"Eᴄ = {self.e_c_ev / self.KB_EV_K:.4f} K·kᴮ = "
              f"{self.e_c_hz / 1e9:.4f} GHz")
        print(f"Eⱼ/Eᴄ = {self.ej_ec_ratio:.2f}")
        print(f"Rₙ = {self.r_n_ohm / 1e3:.1f} kΩ")
        print(f"Δ = {self.delta_material / self.KB_EV_K * 1e3:.3f} mK·kᴮ")
        print(f"Resonator frequency = {resonator_freq_hz / 1e9:.2f} GHz")
        print(f"Coupling g = {coupling_g_hz / 1e6:.1f} MHz")
        print(f"FWHM = {resonator_freq_hz / 10000 / 1e6:.2f} MHz")
        print("=" * 60)
        
        figs = []
        
        # Figure 1: Energy levels
        print("\n[1/4] Plotting energy levels...")
        fig1, _ = self.plot_energy_levels(
            offset_charges, num_levels, resonator_freq_hz, coupling_g_hz
        )
        figs.append(fig1)
        
        # Figure 2: Matrix elements
        print("[2/4] Plotting matrix elements...")
        fig2, _ = self.plot_matrix_elements(
            None, coupling_g_hz, resonator_freq_hz, num_levels
        )
        figs.append(fig2)
        
        # Figure 3: Dispersive shift
        print("[3/4] Plotting dispersive shift...")
        fig3, _ = self.plot_dispersive_shift(
            None, coupling_g_hz, resonator_freq_hz, num_levels
        )
        figs.append(fig3)
        
        # Figure 4: Parity shift vs frequency
        print("[4/4] Plotting parity shift vs frequency...")
        fig4, _ = self.plot_parity_shift_vs_frequency(
            None, coupling_g_hz, num_levels
        )
        figs.append(fig4)
        
        print("\nAll plots complete!")
        return figs


def main():
    """Example usage of OCS class"""
    # Example from the MATLAB script (WashU parameters)
    ej_ec_ratio = 12
    e_j_hz = 8.335e9  # ~8.3 GHz (0.4 K·kB)
    e_c_hz = e_j_hz / ej_ec_ratio  # Hz
    
    # Create OCS instance
    ocs = OCS(
        e_j_hz=e_j_hz,
        e_c_hz=e_c_hz,
        temperature_k=0.012,  # 12 mK
        r_n_ohm=27e3  # 27 kΩ
    )
    
    # Generate all plots
    figs = ocs.plot_all(
        coupling_g_hz=150e6,  # 150 MHz
        resonator_freq_hz=7.0e9,  # 7 GHz
        num_levels=5
    )
    
    plt.show()
    
    return ocs, figs


if __name__ == "__main__":
    ocs, figs = main()

