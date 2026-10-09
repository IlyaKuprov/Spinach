# kernel/states/unit_state.m

- Signature: `rho=unit_state(spin_system)`

## Purpose

Constructs the unit-state representation for the basis already stored in `spin_system`. Call `basis` first; the function uses its formalism and basis data.

## Representations

- In `sphten-liouv`, returns a sparse direct-sum column vector with every substance unit coordinate `bas.offsets(n)+1` set to `chem.concs(n)`, including zero concentrations. These coordinates represent `T(0,0)`.
- In `zeeman-liouv`, vectorises the identity on the spin Hilbert space using MATLAB column-major `(:)` ordering and divides by its Euclidean 2-norm, then weights it by the local concentration and stacks the independent blocks.
- In `zeeman-hilb`, returns a block diagonal matrix of local sparse identities multiplied by their concentrations.
- Other formalism values raise an error.

Thus, “unit” does not mean a trace-one identity in every representation: the two Zeeman cases differ in both shape and normalisation.

## Inputs and output

- `spin_system`: Spinach data object with basis information.
- `rho`: the vector or matrix representation selected above.

## Links

- Basis setup: [`basis.m` page](../basis.md).
- Source: [`kernel/states/unit_state.m`](https://github.com/IlyaKuprov/Spinach/blob/main/kernel/states/unit_state.m).
- Wiki: [`unit_state.m`](https://spindynamics.org/wiki/index.php?title=unit_state.m).

Wavefunction concentration-weighted units are explicitly rejected by `Spinach:unit_state:wavefunction`.

Absent basis metadata or a missing formalism field is rejected by the existing explicit input-validation error before any formalism-specific capability check.
