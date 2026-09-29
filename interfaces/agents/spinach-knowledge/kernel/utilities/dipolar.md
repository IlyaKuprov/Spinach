# kernel/utilities/dipolar.m

**Source:** [https://github.com/IlyaKuprov/Spinach/blob/main/kernel/utilities/dipolar.m](https://github.com/IlyaKuprov/Spinach/blob/main/kernel/utilities/dipolar.m)

## Purpose

Computes dipolar couplings in the presence or absence of periodic boundary conditions. This is an auxiliary function of the Spinach kernel; direct calls are discouraged. The header directs users to `xyz2dd` and `xyz2hfc` to convert Cartesian coordinates into dipolar and hyperfine couplings, respectively.

## Behaviour

- Calls `grumble(spin_system)` to verify that the fields `comp`, `inter`, `chem`, and `tols` are present, erroring with `'spin_system object is missing essential information.'` otherwise.
- Reports the dipolar interaction distance threshold (`spin_system.tols.prox_cutoff`, in Angstrom) and that a dipolar interaction network analysis is running.
- Preallocates an `nspins`-by-`nspins` cell array of distance vectors and loops over chemical species (`spin_system.chem.parts`), then over all ordered pairs of spins within each species.
- For each pair `(n, k)` with both coordinates specified and `n ~= k`, distance vectors from `n` to `k` are determined depending on `spin_system.inter.pbc`:
  - Empty PBC: a single vector `coordinates{k} - coordinates{n}`.
  - Scalar (1D) PBC: `(2*dd_ncells+1)` vectors, one per integer translation `p` in `-dd_ncells:dd_ncells` applied along `pbc{1}`.
  - 2D PBC: `(2*dd_ncells+1)^2` vectors over translations `p` and `q` along `pbc{1}` and `pbc{2}`.
  - 3D PBC: `(2*dd_ncells+1)^3` vectors over translations `p`, `q`, and `r` along `pbc{1}`, `pbc{2}`, and `pbc{3}`.
  - Any other PBC dimension errors with `'PBC translation vector array has invalid dimensions.'`.
- Distance vectors longer than or equal to `spin_system.tols.prox_cutoff` are discarded. If any remaining vector has norm below 0.5 Angstrom, the function errors with a collision message naming spins `n` and `k` (or their PBC images). Surviving vectors are stored in the cell array.
- Builds a sparse proximity matrix `spin_system.inter.proxmatrix` marking non-empty distance-vector cells, finds its nonzero entries, and reports the number of spin pairs under the threshold as `numel(rows)/2` (the matrix is symmetric, so each pair is counted twice).
- For each interacting pair and each PBC distance vector:
  - Computes the distance and the unit vector (`ort`) along the distance vector.
  - Computes the prefactor `A = 0.5 * gamma_n * gamma_k * hbar * mu0 / (4*pi*(distance*1e-10)^3)`, where the 0.5 factor compensates for double counting.
  - Forms the 3-by-3 dipolar coupling matrix `D = A * (I - 3*ort*ort')` written elementwise in the source.
  - If `'sodd'` is listed in `spin_system.sys.enable`, applies an approximate spin-orbit correction by sandwiching `D` between `spin_system.inter.zeeman.ddscal{rows(n)}` and `spin_system.inter.zeeman.ddscal{cols(n)}`.
  - Removes the isotropic part via `D = D - eye(3)*trace(D)/3` to clean up numerical noise.
  - Accumulates `D` into `spin_system.inter.coupling.matrix{rows(n), cols(n)}`, initialising the cell if empty and summing otherwise.

## Inputs and outputs

**Syntax:** `spin_system = dipolar(spin_system)`

**Input:**

- `spin_system` — Spinach data object containing information about chemical subsystems, atomic coordinates, and periodic boundary conditions.

**Output:**

- `spin_system` — Spinach data object with the interaction arrays updated with dipolar and hyperfine coupling information.

## References

- Spinach Wiki page for `dipolar.m`: [https://spindynamics.org/wiki/index.php?title=dipolar.m](https://spindynamics.org/wiki/index.php?title=dipolar.m)
