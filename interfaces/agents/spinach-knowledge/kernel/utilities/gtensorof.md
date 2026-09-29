# kernel/utilities/gtensorof.m

## Purpose

Returns the g-tensor of a specified spin at the input orientation, as documented in the file header. The function is part of the Spinach kernel utilities.

Source: <https://github.com/IlyaKuprov/Spinach/blob/main/kernel/utilities/gtensorof.m>

## Behaviour

- Syntax: `g=gtensorof(spin_system,spin_number)`.
- The function first calls an internal consistency-checking subfunction `grumble(spin_system,spin_number)`.
- The g-tensor is computed as:

  `g = -spin_system.inter.zeeman.ddscal{spin_number} * spin_system.inter.gammas(spin_number) * spin_system.tols.hbar / spin_system.tols.muB`

- The header notes that the same convention (`mu = -mu_b*g*S/hbar`) is used for nuclei, meaning their g-tensors are much smaller than those of electrons.
- Consistency checks performed by `grumble`:
  - Errors with `'spin system object does not contain the required information.'` if `spin_system` does not contain both `inter` and `tols` fields.
  - Errors with `'spin_number must be a positive real integer.'` if `spin_number` is not numeric, not real, not a scalar, is less than 1, or is not an integer.
  - Errors with `'spin_number exceeds the number of spins in the system.'` if `spin_number` is greater than `spin_system.comp.nspins`.

## Inputs and outputs

Inputs:

- `spin_system` — the spin system object; must contain the fields `inter` and `tols` (with `inter.zeeman.ddscal`, `inter.gammas`, `tols.hbar`, `tols.muB`, and `comp.nspins` used in the computation or checks).
- `spin_number` — a positive integer specifying the number of the spin in the `sys.isotopes` list; must be a positive real integer scalar not exceeding the number of spins in the system.

Outputs:

- `g` — a dimensionless 3×3 g tensor; the implementation divides the magnetic-moment coupling by the Bohr magneton, which enters separately in `mu = -mu_B*g*S/hbar`.

## References

- Spinach Wiki page for `gtensorof.m`: <https://spindynamics.org/wiki/index.php?title=gtensorof.m>
- Source file: <https://github.com/IlyaKuprov/Spinach/blob/main/kernel/utilities/gtensorof.m>
