# kernel/reduce.m

- MATLAB source: https://github.com/IlyaKuprov/Spinach/blob/main/kernel/reduce.m
- Wiki: https://spindynamics.org/wiki/index.php?title=reduce.m
- Signature: `projectors=reduce(spin_system,L,rho)`

## Purpose

Returns projectors into independently evolving reduced subspaces, selected using the supplied Liouvillian `L` and state `rho`. This is state-space reduction; this function does not construct a relaxation generator or define relaxation-term units.

## Reduction path

The source first checks whether trajectory-level reduction is disabled by `spin_system.sys.disable` containing `'trajlevel'`; if so, it reports the setting and returns the unit projector `1`. Otherwise the available operations depend on `spin_system.bas.formalism` and the disable settings.

For `zeeman-hilb` and `zeeman-wavef`, the code uses supplied permutation-symmetry irreducible-representation projectors when symmetry treatment is enabled and that information is available. Zero-dimensional irreps are dropped; the state contribution is also screened against `spin_system.tols.irrep_drop`. These formalisms use symmetry screening rather than the Liouville-space zero-track and path-tracing stages.

For `zeeman-liouv` and `sphten-liouv`, the code tries symmetry factorisation when available and not disabled, then applies zero-track elimination when enabled and path tracing to identify disconnected subspaces. Disabling symmetry skips that factorisation; the zero-track and path-tracing stages remain part of this formalism's route, but zero-track elimination returns an identity projector unless `'zte'` is present in `sys.enable`. The detailed reductions therefore depend on the chosen formalism, supplied symmetry data, input state, and configured tolerances.

## Inputs and returned projectors

- `L` - Liouvillian matrix
- `rho` - initial state for source-state screening, or destination state for destination-state screening
- `projectors` - cell array of projectors into the selected subspaces

Use each projector `P` as documented by the source: `L_reduced=P'*L*P` for matrices and `rho_reduced=P'*rho` for state vectors.

## References

- http://dx.doi.org/10.1016/j.jmr.2008.08.008
- http://dx.doi.org/10.1063/1.3398146
- http://dx.doi.org/10.1016/j.jmr.2011.03.010
