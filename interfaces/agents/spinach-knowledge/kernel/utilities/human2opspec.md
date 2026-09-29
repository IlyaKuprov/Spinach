# kernel/utilities/human2opspec.m

## Purpose

Converts user-friendly descriptions of spin states and operators into the formal operator specification (`opspec`) used by the Spinach kernel.

Source: <https://github.com/IlyaKuprov/Spinach/blob/main/kernel/utilities/human2opspec.m>

## Behaviour

The function supports three calling conventions:

1. **String operator and string spin selection** — e.g. `human2opspec(spin_system,'Lz','13C')` returns a list of single-spin opspecs for all spins with the specified isotope name. Valid spin labels are standard isotope names, as well as `'electrons'` (any isotope string starting with `E`), `'nuclei'` (any isotope string not starting with `E`), and `'all'` (all spins in the system). If no matching spins exist, the function errors with `no such spins in the system.`

2. **String operator and numeric spin vector** — e.g. `human2opspec(spin_system,'Lz',[1 2 4])` returns a list of single-spin opspecs for the spins with the given numbers.

3. **Cell array of operator strings and cell array of spin numbers** — e.g. `human2opspec(spin_system,{'Lz','Ly'},{1,2})` returns a product operator specification with `Lz` on spin 1 and `Ly` on spin 2.

Valid operator/state labels for spin-type particles are:

- `'E'` — identity (unit operator)
- `'Lz'`, `'Lx'`, `'Ly'` — Cartesian angular momentum projections
- `'L+'`, `'L-'` — raising and lowering operators
- `'Tl,m'` — irreducible spherical tensor, where `l` and `m` are integers
- `'CTx'`, `'CTy'`, `'CTz'`, `'CT+'`, `'CT-'` — central transition operators in the Zeeman basis
- `'ZLn'` — specific Zeeman energy level projector, where `n` is a positive integer level number

For cavity, lattice, and transmon particle types (`'C'`, `'V'`, `'T'`), the valid labels are `'E'` (unit operator), `'BLn'` (bosonic energy level projector with positive integer level number), and other operator strings passed to `bos2ist` for spherical tensor expansion.

For spin-type particles, the opspec encoding uses index 0 for the unit operator, 1 for the raising operator, 2 for the z projection, and 3 for the lowering operator. The `L+` operator carries a coefficient of `-sqrt(2)` (from the T(1,+1) relationship), `L-` carries `+sqrt(2)`, `Lx` expands into two terms with coefficients `[-sqrt(2); sqrt(2)]/2`, and `Ly` expands into two terms with coefficients `[-sqrt(2); -sqrt(2)]/2i`. Irreducible spherical tensors are encoded via `lm2lin(l,m)` with validation that `l >= 0`, `|m| <= l`, both are integers, and `lm2lin(l,m)+1` does not exceed the square of the spin multiplicity. Zeeman level projectors are expanded via `enlev2ist` with the `'S'` particle flag; bosonic level projectors use the `'B'` flag. Central transition operators are expanded via `ct2ist` using the spin multiplicity.

Input validation (in the internal `grumble` function) enforces that operator/spins argument type combinations are valid, that cell arrays of operators and spins have equal element counts, that all operator cell elements are strings, that numeric spin lists are real positive integer rows without repetitions, that cell-array spin entries are positive integers without repetitions, and that the spin list is not empty.

The function notes that direct calls are not necessary; users should call `operator.m` and `state.m` instead.

## Inputs and outputs

**Inputs:**

- `spin_system` — Spinach spin system structure
- `operators` — operator/state label string, or cell array of label strings
- `spins` — isotope name or group string (`'all'`, `'electrons'`, `'nuclei'`), numeric row vector of spin numbers, or cell array of spin numbers

**Outputs:**

- `opspecs` — Spinach operator specification: a cell array of row vectors specifying which operator enters the Kronecker product for which spin
- `coeffs` — coefficient with which each of the Kronecker products enters the overall sum

## References

- Spinach Wiki page for this function: <https://spindynamics.org/wiki/index.php?title=human2opspec.m>
- Source code: <https://github.com/IlyaKuprov/Spinach/blob/main/kernel/utilities/human2opspec.m>
