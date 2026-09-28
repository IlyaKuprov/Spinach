# kernel/utilities/path_trace.m

- Signature: `projectors=path_trace(spin_system,L,rho)`

## Purpose

`path_trace` treats the supplied Liouvillian as a graph adjacency matrix, finds its weakly connected subgraphs, and returns projectors into independently evolving populated subspaces.

## Parameters

- `spin_system` — supplies run settings, tolerances, and basis formalism.
- `L` — Hamiltonian or Liouvillian matrix.
- `rho` — initial state for source-state screening or detection state for destination-state screening. Pass `[]` to disable population screening.

## Outputs

`projectors` is a cell array of projectors. For a projector `P`, use `L_reduced=P'*L*P` for matrices and `rho_reduced=P'*rho` for state vectors.

## Algorithm and run conditions

The function requires numeric `L` and `rho`, a square `L`, and, when `rho` is nonempty, matching dimensions between the columns of `L` and the rows of `rho`. If `pt` is in `spin_system.sys.disable`, or if `size(L,2)<spin_system.tols.merge_dim`, it skips path tracing and returns the unit projector `{1}`.

Otherwise, it forms a connectivity matrix from `abs(L)>spin_system.tols.liouv_zero`, symmetrizes it, adds the identity to retain isolated states, and finds connected components with `scomponents`. All components are retained when `rho` is empty. With a nonempty `rho`, components are retained when their population exceeds `spin_system.tols.subs_drop`: `sphten-liouv` and `zeeman-liouv` use the 1-norm of the corresponding entries of the state vector; `zeeman-hilb` checks both corresponding rows and columns of the state matrix. An unexpected formalism raises an error.

The function builds sparse projectors for retained components and reports their dimensions. Unless `merge` is in `spin_system.sys.disable`, it groups small subspaces using `binpack(subspace_dims,spin_system.tols.merge_dim)` and concatenates their projectors into working subspaces.

## Further information

- http://dx.doi.org/10.1063/1.3398146
- http://dx.doi.org/10.1016/j.jmr.2011.03.010
- https://spindynamics.org/wiki/index.php?title=path_trace.m