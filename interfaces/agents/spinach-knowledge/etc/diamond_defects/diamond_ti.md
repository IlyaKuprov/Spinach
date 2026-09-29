# etc/diamond_defects/diamond_ti.m

- MATLAB implementation: [etc/diamond_defects/diamond_ti.m](https://github.com/IlyaKuprov/Spinach/blob/main/etc/diamond_defects/diamond_ti.m)

**Call:** `[sys,inter] = diamond_ti(parameters)`

Builds the N3 or OK1 titanium-related defect model in diamond, with a nitrogen spin, optional titanium isotope and (for OK1 only) up to two nearby `13C` spins.

## Inputs

- `parameters.centre`: `'n3'` or `'ok1'`; matching is case-insensitive.
- `parameters.orientation`: `'111'`, `'110'` or `'100'`; the selected crystal direction is aligned with the magnetic field. These values are matched exactly.
- `parameters.titanium`: isotope label, or `'none'` to omit titanium.
- `parameters.n_13c`: for OK1, use a count of 0, 1 or 2. The source checks its range but does not test integrality. N3 is set to zero internally; a supplied nonzero count is rejected for N3.

The function uses degrees for its internal frame tilts (`cosd`/`sind` rotations).

## Centre-specific parameters

The source assigns these principal `g` values and nitrogen (`An`) and titanium (`Ati`) hyperfine values:

| Centre | `g` values | `An` (mT) | `Ati` (mT) | `g`-frame tilt | hyperfine-frame tilt |
| --- | --- | --- | --- | ---: | ---: |
| N3 | [2.0022, 2.0025, 2.0020] | [0.11, 0.15, 0.11] | [0.28, 0.40, 0.28] | 32° | 26° |
| OK1 | [2.0031, 2.0019, 2.0025] | [0.55, 0.77, 0.54] | [0.06, 0.06, 0.06] | 40° | 20° |

The nitrogen isotope is fixed to `14N`. If titanium is not `'none'`, the same centre-specific `Ati` tensor is assigned to the requested isotope. OK1 carbons use principal hyperfine values [2.62, 2.62, 4.38] mT and two specified local frames; `n_13c` selects the first zero, one or two of those frames. Conversion from mT uses `abs(spin('E'))/(2*pi)*1e-3`.

## Outputs and reference

- `sys`: electron, `14N`, optional titanium isotope and requested OK1 carbons.
- `inter`: electron Zeeman tensor and anisotropic electron–nuclear coupling matrices.

Magnetic parameters are attributed to Nadolinny et al., *Crystals* **7**, 237 (2017), [doi:10.3390/cryst7080237](https://doi.org/10.3390/cryst7080237). The routine is not a coordinate-relaxation model.

[Spin Dynamics Wiki source page](https://spindynamics.org/wiki/index.php?title=diamond_ti.m).
