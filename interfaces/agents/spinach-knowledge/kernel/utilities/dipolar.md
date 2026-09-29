# kernel/utilities/dipolar.m

**Source:** [https://github.com/IlyaKuprov/Spinach/blob/main/kernel/utilities/dipolar.m](https://github.com/IlyaKuprov/Spinach/blob/main/kernel/utilities/dipolar.m)

## Purpose

Computes dipolar couplings in the presence or absence of periodic boundary conditions. This is an auxiliary function of the Spinach kernel; direct calls are discouraged. The header directs users to `xyz2dd` and `xyz2hfc` to convert Cartesian coordinates into dipolar and hyperfine couplings, respectively.

## Physical construction

For coordinate-bearing spin pairs within each chemical species, an isolated system contributes one separation vector. With periodic boundary conditions, integer-translated images contribute along one, two or three lattice vectors, from `-dd_ncells` to `+dd_ncells` in each periodic direction. Only separations shorter than `spin_system.tols.prox_cutoff` (ångström) enter the proximity network; a separation below `0.5` Å is treated as an atomic collision.

Each retained image contributes a dipolar tensor proportional to `gamma_n*gamma_k/r^3` and to the traceless axial form given by the identity matrix minus three times the outer product of the unit separation vector with itself; `r` is converted from ångström to metres. The prefactor contains `0.5` to compensate for counting ordered spin pairs twice. When the `sodd` spin-orbit option is enabled, Zeeman anisotropy scaling acts on the tensor; its numerical isotropic trace is then removed. Contributions from accepted images sum into the pair’s dipolar coupling matrix.

## Inputs and outputs

**Syntax:** `spin_system = dipolar(spin_system)`

**Input:**

- `spin_system` — Spinach data object containing information about chemical subsystems, atomic coordinates, and periodic boundary conditions.

**Output:**

- `spin_system` — Spinach data object with the interaction arrays updated with dipolar and hyperfine coupling information.

## References

- Spinach Wiki page for `dipolar.m`: [https://spindynamics.org/wiki/index.php?title=dipolar.m](https://spindynamics.org/wiki/index.php?title=dipolar.m)
