# Material System Quick Reference

## Installation

```bash
pip install pyyaml>=5.4.0  # Required for materials.yaml
```

## Quick Start

```python
from ocs_transmon import OCS

# Default (aluminum)
ocs = OCS(8.335e9, 0.695e9)

# Specify material
ocs = OCS(8.335e9, 0.695e9, material='hafnium')

# List available
print(OCS.list_materials())
```

## Available Materials

| Material | Tc [K] | Δ [μeV] | Common Use |
|----------|--------|---------|------------|
| aluminum | 1.200 | 189 | Standard transmons |
| hafnium | 0.128 | 22.5 | Low-temp studies |
| aluminum_manganese | 0.100 | 15.7 | Very low Tc |
| niobium | 9.200 | 1400 | Resonators |
| titanium_nitride | 3.000 | 440 | High-Q devices |

## Common Operations

```python
# List materials
OCS.list_materials()

# Get properties
props = OCS.get_material_properties('hafnium')

# Override properties
ocs = OCS(8e9, 0.7e9, material='aluminum', tc=0.8)

# Custom gaps
ocs = OCS(8e9, 0.7e9, delta_l_hz=50e9, delta_r_hz=45e9)

# Access material info
print(f"{ocs.material_name}: Tc={ocs.tc}K, DOS={ocs.dos:.2e}")
```

## Add New Material

Edit `materials.yaml`:

```yaml
materials:
  new_material:
    name: "Material Name"
    symbol: "Sym"
    description: "Description"
    properties:
      dos: 2.0e10
      fermi_level: 10.0
      tc: 1.5
      delta: 2.2e-4
      notes: "Notes here"
```

## Examples

```bash
python example_materials.py  # Material comparison
python verify_translation.py  # Verify still works
```

## Breaking Changes

Removed class constants (use instance variables instead):
- ~~`OCS.DOS_AL`~~ → `ocs.dos`
- ~~`OCS.TC_AL`~~ → `ocs.tc`
- ~~`OCS.DELTA_AL`~~ → `ocs.delta_material`

## Files

- `materials.yaml` - Material database
- `MATERIALS_GUIDE.md` - Full documentation
- `example_materials.py` - Usage examples
- `MATERIALS_UPDATE_SUMMARY.md` - Implementation details

