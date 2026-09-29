# kernel/utilities/cubic_lattice.m

## Purpose

Creates a periodic volume-centred cubic lattice with user-supplied parameters, returning Spinach input data structures suitable for periodic-boundary simulations.

Source: <https://github.com/IlyaKuprov/Spinach/blob/main/kernel/utilities/cubic_lattice.m>

## Behaviour

- Syntax: `[sys,inter]=cubic_lattice(isotope,spacing,n_periods)`.
- Calls the internal consistency-checking function `grumble` on the inputs before building the lattice.
- Builds the isotope list as a `1 x n_periods^3` cell array, with every element set to the supplied isotope string.
- Builds `inter.coordinates` as an `n_periods^3 x 1` cell array of coordinate triplets. For indices `n`, `k`, `m` each running from 1 to `n_periods`, the coordinate at linear index `sub2ind([n_periods n_periods n_periods],m,k,n)` is set to `spacing*[(n-1) (k-1) (m-1)]`.
- Builds periodic boundary translation vectors in `inter.pbc` as `spacing*n_periods*[1 0 0]`, `spacing*n_periods*[0 1 0]`, and `spacing*n_periods*[0 0 1]`.
- Input validation (in `grumble`):
  - `isotope` must be a character string, otherwise the function errors with `isotope must be a character string.`.
  - `spacing` must be a scalar, real, finite, positive number, otherwise the function errors with `spacing must be a positive real number.`.
  - `n_periods` must be a scalar, real, finite number that is at least 1 and an integer (i.e. `mod(n_periods,1)==0`), otherwise the function errors with `n_periods must be a positive real integer.`.

## Inputs and outputs

Inputs:

- `isotope` — character string specifying the isotope, for example `'13C'`.
- `spacing` — lattice spacing in Angstroms.
- `n_periods` — number of lattice periods in each of the three spatial dimensions.

Outputs:

- `sys`, `inter` — Spinach input data structures with the following fields set: `sys.isotopes`, `inter.coordinates`, `inter.pbc`.

## References

- Spinach Wiki: <https://spindynamics.org/wiki/index.php?title=cubic_lattice.m>
- Source file: <https://github.com/IlyaKuprov/Spinach/blob/main/kernel/utilities/cubic_lattice.m>
