# etc/diamond_defects/diamond_ti.m

`[sys,inter]=diamond_ti(parameters)`

Builds Spinach models for the N3 and OK1 titanium-related centres. The magnetic parameters are attributed to Nadolinny et al., *Crystals* **7**, 237 (2017) ([doi:10.3390/cryst7080237](https://doi.org/10.3390/cryst7080237)).

## Parameters

- `parameters.centre`: `'n3'` or `'ok1'` (case-insensitive).
- `parameters.orientation`: `'111'`, `'110'`, or `'100'`; the specified crystal plane normal is aligned with the field axis.
- `parameters.titanium`: titanium isotope label, or `'none'` to omit titanium.
- `parameters.n_13c`: for OK1, an integer from 0 to 2 is required and selects the reported nearest-neighbour carbon-13 sites. It is not supported for N3; omit it or set it to zero.

Both models include an electron labelled `E` and a nitrogen-14 nucleus. Principal g values and hyperfine values (in mT, converted internally to frequency units) are:

| Centre | g | ¹⁴N hyperfine | Ti hyperfine |
|---|---|---|---|
| N3 | 2.0022, 2.0025, 2.0020 | 0.11, 0.15, 0.11 | 0.28, 0.40, 0.28 |
| OK1 | 2.0031, 2.0019, 2.0025 | 0.55, 0.77, 0.54 | 0.06, 0.06, 0.06 |

The titanium interaction is included only when an isotope other than `'none'` is requested. For OK1, each selected carbon-13 site has principal hyperfine values 2.62, 2.62, and 4.38 mT. The routine rotates the centre-specific tensor frames into the requested field orientation.

The centre, orientation, and titanium label must be character values; the orientation must be one of the three listed values. For OK1, `n_13c` is required and limited to 0–2; a nonzero value is rejected for N3.

[Spin Dynamics Wiki page](https://spindynamics.org/wiki/index.php?title=diamond_ti.m).
