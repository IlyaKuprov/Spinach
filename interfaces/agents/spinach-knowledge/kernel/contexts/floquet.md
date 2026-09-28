# kernel/contexts/floquet.m

- Signature: `[answer,sph_grid]=floquet(spin_system,pulse_sequence,...`

## Purpose

Floquet magic-angle-spinning context. It builds a Liouvillian superoperator for a powder simulation and passes it to a pulse sequence supplied as a function handle. The full call is `[answer,sph_grid]=floquet(spin_system,pulse_sequence,parameters,assumptions)`. The `assumptions` character string is passed to `assume` when the Hamiltonian is built.

## Physical / mathematical content

The context represents rotor-periodic dynamics in Floquet space. For each orientation in a spherical averaging grid, it constructs Fourier terms from the Hamiltonian's non-empty spherical interaction ranks, adds the isotropic Hamiltonian and the rotor turning generator, and calls the pulse sequence with the resulting Floquet Liouvillian, relaxation operator, and kinetics operator.

The Floquet spatial dimension is `2*parameters.max_rank+1`; the spin dimension is the dimension of the isotropic Hamiltonian. The retained harmonic rank should be increased until the result converges. Slower spinning generally requires more ranks, and the required rank is approximately the number of spinning sidebands. The code warns if `max_rank` is below a non-empty interaction rank because those harmonics are truncated.

## Parameters and constraints

- `parameters.rate`: spinning rate in Hz. Positive values correspond to the JEOL spinning direction; negative values correspond to the Varian and Bruker direction. Its spinning sense matches `singlerot.m`: the same rate produces the same powder result in both contexts.
- `parameters.axis`: spinning axis as a normalized three-element vector. The implementation requires a row vector of three real numbers.
- `parameters.spins`: non-empty cell array of spin isotope names involved in the pulse sequence, such as `{'1H','13C'}`. Each isotope must occur in the spin system.
- `parameters.offset`: transmitter offsets in Hz corresponding to `parameters.spins`. If omitted, zero offsets are used.
- `parameters.max_rank`: maximum harmonic rank retained in the Floquet calculation; required to be a non-negative integer.
- `parameters.grid`: filename of a spherical grid in `kernel/grids`. Single-crystal simulations are not supported; use `singlerot.m` instead.
- `parameters.sum_up`: when `1` (default), returns the weighted powder average; when `0`, returns the pulse-sequence result for each orientation as a cell array.
- Pulse sequences may require additional fields in `parameters`; consult the pulse sequence documentation. The context also supplies `parameters.spc_dim` and `parameters.spn_dim`, the spatial and spin dynamics subspace dimensions, to the pulse sequence.

The context requires the `zeeman-liouv` or `sphten-liouv` Liouville-space formalism and does not support a numerical rotating-frame transformation through `parameters.rframes`. Perturbative corrections to the rotating-frame transformation are not supported; use `singlerot.m` instead.

## Numerical / algorithmic content

The context loads grid angles and weights, evaluates each powder orientation, and forms a weighted sum unless `parameters.sum_up` is disabled. It supports parallel orientation evaluation through MATLAB's parallel computing facilities; `parameters.serial` can turn parallel execution off. The state projector assumes a powder, not a single crystal.

## Outputs

- `answer`: the weighted powder average of the pulse-sequence output, or a cell array of individual orientation outputs when `parameters.sum_up` is `0`.
- `sph_grid`: the spherical grid used in the calculation, containing `alphas`, `betas`, `gammas`, and `weights`.

<https://spindynamics.org/wiki/index.php?title=floquet.m>