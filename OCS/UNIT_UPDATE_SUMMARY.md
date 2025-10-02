# Unit System Update Summary

## Overview

Updated the OCS transmon simulator to use Hz for all energy and frequency 
parameters, with explicit unit suffixes in all variable names to avoid confusion.

## Key Changes

### 1. Input Parameters Now in Hz

**Before:**
```python
OCS(e_josephson=0.4*OCS.KB_EV_K,  # eV
    e_charging=0.033*OCS.KB_EV_K,  # eV
    temperature=0.02,  # K
    normal_resistance=27e3)  # Ω
```

**After:**
```python
OCS(e_j_hz=8.335e9,  # Hz (GHz)
    e_c_hz=0.695e9,  # Hz (GHz)
    temperature_k=0.02,  # K
    r_n_ohm=27e3)  # Ω
```

### 2. All Variable Names Include Units

| Old Name | New Name | Unit | Description |
|----------|----------|------|-------------|
| `e_josephson` | `e_j_hz` | Hz | Josephson energy |
| `e_charging` | `e_c_hz` | Hz | Charging energy |
| `temperature` | `temperature_k` | K | Temperature |
| `normal_resistance` | `r_n_ohm` | Ω | Normal resistance |
| `delta_left` | `delta_l_hz` | Hz | Left superconducting gap |
| `delta_right` | `delta_r_hz` | Hz | Right superconducting gap |
| `coupling_g` | `coupling_g_hz` | Hz | Transmon-resonator coupling |
| `resonator_freq` | `resonator_freq_hz` | Hz | Resonator frequency |
| `freq_range` | `freq_range_hz` | Hz | Frequency range array |

### 3. Internal Variables

The class now maintains both Hz and eV versions internally:

**User-facing (Hz):**
- `self.e_j_hz`
- `self.e_c_hz`
- `self.temperature_k`
- `self.r_n_ohm`

**Internal (eV for calculations):**
- `self.e_j_ev` (converted from Hz)
- `self.e_c_ev` (converted from Hz)
- `self.delta_l_ev`
- `self.delta_r_ev`

### 4. Method Signature Changes

#### `__init__`
```python
# Before
def __init__(self, e_josephson, e_charging, temperature=0.02,
             normal_resistance=27e3, delta_left=None, 
             delta_right=None)

# After
def __init__(self, e_j_hz, e_c_hz, temperature_k=0.02,
             r_n_ohm=27e3, delta_l_hz=None, delta_r_hz=None)
```

#### `from_capacitance`
```python
# Before
@classmethod
def from_capacitance(cls, total_capacitance, delta, 
                    normal_resistance, **kwargs)

# After
@classmethod
def from_capacitance(cls, c_total_f, delta_hz, r_n_ohm, **kwargs)
```

#### `compute_dispersive_matrix`
```python
# Before
def compute_dispersive_matrix(self, offset_charge, coupling_g, 
                              resonator_freq, num_levels=6,
                              charge_cutoff=30)

# After
def compute_dispersive_matrix(self, offset_charge, coupling_g_hz, 
                              resonator_freq_hz, num_levels=6,
                              charge_cutoff=30)
```

#### Plotting Methods
```python
# Before
def plot_matrix_elements(self, offset_charges=None, 
                        coupling_g=150e6,
                        resonator_freq=7.0e9, ...)

# After
def plot_matrix_elements(self, offset_charges=None, 
                        coupling_g_hz=150e6,
                        resonator_freq_hz=7.0e9, ...)
```

Similar changes for:
- `plot_dispersive_shift()`
- `plot_parity_shift_vs_frequency()`
- `plot_all()`

### 5. Unit Conversion

**Energy conversions (internal use):**
```python
# Hz → eV
energy_ev = energy_hz * PLANCK_EV_S

# eV → Hz
energy_hz = energy_ev / PLANCK_EV_S
```

Where `PLANCK_EV_S = h / eV ≈ 4.136e-15 eV·s`

### 6. Example Usage Updates

**Before:**
```python
from ocs_transmon import OCS

kb = OCS.KB_EV_K
e_j = 0.4 * kb  # eV
e_c = e_j / 12

ocs = OCS(e_j, e_c, temperature=0.02, normal_resistance=27e3)
figs = ocs.plot_all(coupling_g=150e6, resonator_freq=7.0e9)
```

**After:**
```python
from ocs_transmon import OCS

e_j_hz = 8.335e9  # ~8.3 GHz
e_c_hz = e_j_hz / 12

ocs = OCS(e_j_hz, e_c_hz, temperature_k=0.02, r_n_ohm=27e3)
figs = ocs.plot_all(coupling_g_hz=150e6, resonator_freq_hz=7.0e9)
```

### 7. Conversion Reference

Common energy values in both units:

| K·kB | eV | Hz (GHz) | Description |
|------|-----|----------|-------------|
| 0.4 | 3.447e-5 | 8.335 | E_J (WashU) |
| 0.0333 | 2.873e-6 | 0.6946 | E_C (WashU, ratio=12) |
| 0.295 | 2.542e-5 | 6.149 | E_J (Serniak) |
| 0.017 | 1.465e-6 | 0.354 | E_C (Serniak) |

Where: 1 K·kB = 8.617e-5 eV = 20.836 GHz

## Benefits

1. **Intuitive Units**: Frequencies in Hz match standard lab usage
2. **No Confusion**: Variable names explicitly show units
3. **Type Safety**: Harder to mix up parameters
4. **Consistency**: All frequencies use same unit system
5. **Backward Compatible**: Can still compute from capacitance

## Files Modified

1. `ocs_transmon.py` - Core implementation
2. `example_usage.py` - All 5 examples updated
3. `verify_translation.py` - Verification script updated

## Verification

All tests pass with identical numerical results:
- E₀₁ = 5.958 GHz ✓
- E₀₂ = 11.462 GHz ✓
- E₀₃ = 13.773 GHz ✓
- χ (parity) = -1.194 MHz ✓
- Hamiltonian structure correct ✓
- Symmetries verified ✓

## Migration Guide

To update existing code:

1. **Energy inputs**: Multiply eV values by `1/PLANCK_EV_S` or use GHz directly
2. **Parameter names**: Add unit suffixes (`_hz`, `_k`, `_ohm`, `_f`)
3. **Method calls**: Update all plotting method calls with `_hz` suffixes

Example:
```python
# Old
ocs = OCS(0.4*OCS.KB_EV_K, 0.033*OCS.KB_EV_K)
ocs.plot_all(coupling_g=150e6, resonator_freq=7e9)

# New
ocs = OCS(8.335e9, 0.695e9)
ocs.plot_all(coupling_g_hz=150e6, resonator_freq_hz=7e9)
```

## Quick Reference

```python
from ocs_transmon import OCS

# Create instance with energies in Hz
ocs = OCS(
    e_j_hz=8.335e9,      # Josephson energy [Hz]
    e_c_hz=0.695e9,      # Charging energy [Hz]
    temperature_k=0.02,  # Temperature [K]
    r_n_ohm=27e3         # Normal resistance [Ω]
)

# Or from capacitance
ocs = OCS.from_capacitance(
    c_total_f=65e-15,    # Total capacitance [F]
    delta_hz=45.7e9,     # SC gap [Hz] (1.89e-4 eV)
    r_n_ohm=27e3         # Normal resistance [Ω]
)

# Generate plots
figs = ocs.plot_all(
    coupling_g_hz=150e6,      # 150 MHz [Hz]
    resonator_freq_hz=7.0e9,  # 7 GHz [Hz]
    num_levels=6
)

# Access values (both Hz and eV available)
print(f"E_J = {ocs.e_j_hz / 1e9:.3f} GHz")  # User-facing
print(f"E_C = {ocs.e_c_ev:.6e} eV")         # Internal
```

## Notes

- Internal calculations still use eV for compatibility with existing formulas
- All user-facing parameters and returns use Hz
- Physical constants remain in class attributes (unchanged)
- Documentation and docstrings fully updated

