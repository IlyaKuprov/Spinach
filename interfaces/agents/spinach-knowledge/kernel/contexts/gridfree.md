# kernel/contexts/gridfree.m

- Signature: `answer=gridfree(spin_system,pulse_sequence,parameters,assumptions)`

## Purpose

Fokker-Planck magic angle spinning and SLE context. Generates a Liouvillian superoperator and passes it to the pulse sequence function, which should be supplied as a handle.

## Physical / mathematical content

- The file uses a Fokker-Planck-style enlarged state space in which spatial or orientational coordinates are promoted to extra dimensions and coupled to spin dynamics through differential operators.
- The context calls `relaxation(spin_system)` and `kinetics(spin_system)` and includes their operators in the enlarged state space. Rotational correlation times for SLE are specified by `parameters.tau_c`; `inter.tau_c` is used by the Redfield theory module.

## Numerical / algorithmic content

- The context obtains spatial operators from `sle_operators`, combines isotropic and anisotropic spin terms with spatial operators, and adds rotation and rotational diffusion operators when their parameters are specified.

## Parameters / inputs

- `pulse_sequence` — a function handle to a pulse sequence in the experiments directory.
- `assumptions` — a string passed to `assume.m` when the Hamiltonian is built.
- `parameters` — a structure with the following subfields:
  - `.rate` — spinning rate in Hz: positive for JEOL, negative for Varian and Bruker, due to different rotation directions.
  - `.axis` — spinning axis, given as a normalized 3-element vector.
  - `.spins` — a cell array of spins involved in the pulse sequence, e.g. `{'1H','13C'}`.
  - `.offset` — a cell array of transmitter offsets in Hz for the spins listed in `parameters.spins`.
  - `.max_rank` — maximum D-function rank retained in the solution. Increase until convergence is achieved; it is approximately equal to the number of spinning sidebands in the spectrum.
  - `.tau_c` — rotational diffusion correlation times in seconds: a single number for isotropic diffusion or a symmetric positive-definite 3x3 correlation-time tensor for anisotropic diffusion. The rotational diffusion tensor is `inv(6*tau_c)`.
  - `.add_terms` — optional cell array of two-element cell arrays `{c,A}`; each term `c*A` is added to the isotropic spin Hamiltonian after frequency offsets have been applied.
  - `.*` — additional subfields may be required by the pulse sequence; check its documentation.
- The parameters structure is passed to the pulse sequence with `parameters.spc_dim` (matrix dimension of the spatial dynamics subspace) and `parameters.spn_dim` (matrix dimension of the spin dynamics subspace) set.

## Outputs

- Returns the powder average of whatever the pulse sequence returns.
- The Wigner D-function rank truncation level depends on the spinning rate: slower spinning requires greater ranks.
- Rotational correlation times for SLE go in `parameters.tau_c`, not `inter.tau_c` (the latter is used by the Redfield theory module).
- The state projector assumes a powder; single-crystal MAS is not currently supported.
- Perturbative corrections to the rotating-frame transformation are not supported; use `singlerot.m` if needed.

## Implementation structure

The context builds the Liouvillian, projects the initial and detection states into `D[0,0,0]`, and calls `pulse_sequence(spin_system,parameters,H,R,K)`.

<https://spindynamics.org/wiki/index.php?title=gridfree.m>