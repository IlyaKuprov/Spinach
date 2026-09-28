# kernel/reduce.m

- Signature: `projectors=reduce(spin_system,L,rho)`

## Purpose

`reduce` returns projectors onto independently evolving subspaces selected for the supplied state `rho` and Liouvillian `L`. The available reductions depend on the basis formalism and the methods disabled in `spin_system.sys.disable`. For Liouville-space formalisms, the sequence is permutation-symmetry factorisation (when available), zero-track elimination, and disconnected-subspace identification by path tracing.

## Physical / mathematical content

For a projector `P`, the reduced Liouvillian and state vector are `L_reduced=P'*L*P` and `rho_reduced=P'*rho`, respectively.

## Numerical / algorithmic content

Permutation-symmetry projectors with zero dimension or a contribution below `spin_system.tols.irrep_drop` are discarded. In `zeeman-liouv` and `sphten-liouv`, zero-track elimination and path tracing further split the retained subspaces. In `zeeman-hilb` and `zeeman-wavef`, the implementation uses symmetry screening only. If trajectory-level reduction is disabled with `'trajlevel'`, the function returns the unit projector `1`.

## Parameters / inputs

- `L` - Liouvillian matrix.
- `rho` - initial state for source-state screening, or destination state for destination-state screening.

## Outputs

- `projectors` - cell array of projectors into independently evolving reduced subspaces. For a projector `P`, use `L_reduced=P'*L*P` for matrices and `rho_reduced=P'*rho` for state vectors.

## References

- <http://dx.doi.org/10.1016/j.jmr.2008.08.008>
- <http://dx.doi.org/10.1063/1.3398146>
- <http://dx.doi.org/10.1016/j.jmr.2011.03.010>
- <https://spindynamics.org/wiki/index.php?title=reduce.m>
