# OpenMC Reactor Physics Studies

A progressive set of Monte Carlo neutronics calculations built with
[OpenMC](https://openmc.org), exploring the core physics of light-water
reactor lattices. Built as a self-directed learning project to demonstrate
reactor-physics reasoning and Monte Carlo modelling skills.

**Nuclear data:** ENDF/B-VIII.0 (continuous-energy HDF5 library)

## Overview

Each script is a self-contained study, ordered from a bare critical mass
up to a full PWR-type assembly with reactivity-coefficient analysis.

| File | Study | Key result |
|------|-------|-----------|
| `01_sphere.py` | Bare enriched-uranium metal sphere | k-effective vs radius; leakage-dominated criticality |
| `02_cylinder.py` | UO₂ pin cell (fuel / Zr clad / water) | Infinite-lattice k∞ with thermal scattering S(α,β) |
| `03_cylinder_kinf.py` | Moderation study: k∞ vs lattice pitch | Bell-shaped curve; confirms PWR pin is under-moderated |
| `04_cylinder_grid_3x3.py` | 3×3 lattice with a central water hole | Universe nesting; local moderation effect |
| `05_LWR_assembly_17x17.py` | 17×17 PWR-type assembly | Power map, neutron spectrum, void coefficient |

## Physics highlights

### Moderation curve (`03`)

k∞ as a function of pitch shows a maximum near ~1.5 cm. The standard PWR
pitch (1.26 cm) sits to the left of this peak, confirming that LWR lattices
are deliberately *under-moderated* — the basis of a negative moderator/void
coefficient.

![Moderation curve: k-infinity vs pitch](figures/moderation_curve.png)

### Power map (`05`)

A mesh-tallied fission-rate map of the assembly reveals power peaking in the
fuel pins adjacent to the water-filled guide tubes, where local thermalisation
is enhanced. Radial peaking factor ≈ 1.14.

![Assembly power map](figures/power_map_1.00.png)

### Neutron spectrum (`05`)

Flux per unit lethargy shows the three expected regions: the fast fission
bump (~1 MeV), the 1/E slowing-down plateau, and the thermal peak (~0.025 eV).
U-238 resonance dips are visible in the epithermal range.

![Neutron spectrum](figures/spectrum.png)

### Void coefficient (`05`)

Sweeping moderator density from 1.0 to ~0.01 g/cm³ gives a strongly negative
reactivity change, confirming the passive-safety behaviour of under-moderated
light-water lattices. The response is non-linear, reflecting movement along
the moderation curve.

![k-infinity vs moderator density](figures/void_coefficient.png)

## Running

```bash
conda activate openmc-env
python 01_sphere.py
```

`05_LWR_assembly_17x17.py` is controlled by three flags at the top of the file:

- `use_only_cellule` — single pin cell instead of the full assembly
- `calcul_energy_spectrum` — enable spectrum and power-map tallies
- `calculate_density` — sweep moderator density (void-coefficient study)

## Modelling assumptions and limitations

This is a learning project, not a validated production model. Simplifications:

- Cold conditions (water at 1.0 g/cm³, 294 K) — not the hot operating state
- Pure zirconium cladding (not Zircaloy); no pellet–clad gap modelled
- 2D infinite lattice (axial dimension not represented)
- Guide-tube layout based on a standard French PWR assembly
  (ref: *Exploitation des cœurs REP*, Génie Atomique collection, INSTN);
![Geometry LWR 17x17](figures/05_LWR_assembly_17x17_geom.png)
- Results are not benchmarked against reference criticality data

## Requirements

openmc
numpy
matplotlib
Plus an ENDF/B-VIII.0 HDF5 cross-section library with the
`OPENMC_CROSS_SECTIONS` environment variable set.