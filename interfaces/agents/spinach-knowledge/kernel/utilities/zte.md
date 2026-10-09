# kernel/utilities/zte.m

## Purpose

`zte.m` performs zero track elimination: it inspects the first few steps of the system trajectory and drops states that did not get populated beyond a user-specified tolerance, returning a projector matrix into the reduced state space.

## Behaviour

- Syntax: `projector=zte(spin_system,L,rho,nstates)`.
- Input validation (`grumble`) requires the basis formalism to be `zeeman-liouv` or `sphten-liouv`, both `L` and `rho` to be numeric, `rho` to be a single vector (not a stack), `L` to be square, and `size(L,2)==size(rho,1)`.
- If `nstates` is supplied, it must be a real positive integer scalar not exceeding `numel(rho)`; otherwise an error is raised.
- Unless `'zte'` is listed in `spin_system.sys.enable`, the function reports that zero track elimination is not enabled, the basis is left unchanged, and `projector=1` is returned.
- If `nnz(rho)/numel(rho) > spin_system.tols.zte_maxden`, the function skips elimination (too few zeros in the state vector) and returns `projector=1`.
- If `norm(rho,1) < spin_system.tols.zte_tol`, the function skips elimination (state vector norm below drop tolerance, too small for the Krylov procedure) and returns `projector=1`.
- Otherwise, the time step is set to `1/cheap_norm(L)`; if this is infinite (zero Liouvillian), a unit time step is used with a report.
- The trajectory is preallocated as a complex matrix of size `numel(rho)`-by-`spin_system.tols.zte_nsteps`, with `trajectory(:,1)=rho`.
- Steps 2 through `spin_system.tols.zte_nsteps` are computed with the Krylov `step` function. After each step, the active space dimension (number of states whose maximum absolute amplitude over the trajectory exceeds `spin_system.tols.zte_tol`) is compared with the previous value; the loop terminates early when the dimension stops changing.
- Track selection: if `nstates` is given, states are ranked by their maximum absolute amplitude over the trajectory (descending) and the top `nstates` are kept; otherwise all states whose maximum absolute amplitude is below `spin_system.tols.zte_tol` are dropped.
- In the compiled spherical-tensor space, every substance unit coordinate is retained regardless of its trajectory weight. This also applies to zero-population and spin-free substances; these mandatory coordinates may increase the retained dimension beyond `nstates`. For symmetry-reduced calls, `reduce` supplies the support of the projected unit directions instead of the original offsets.
- The projector is built as `speye(size(L))` with the columns corresponding to zero tracks deleted. The intended usage is `L_reduced=P'*L*P` and `rho_reduced=P'*rho`.
- The default tolerance may be altered by setting `sys.tols.zte_tol` before calling `create.m`.
- If tiny interactions or nearly equivalent spins are present, it is best to leave zero track elimination off by omitting `'zte'` from the `sys.enable` cell array.

## Inputs and outputs

Inputs:
- `spin_system` — spin system object supplying tolerances (`zte_tol`, `zte_maxden`, `zte_nsteps`), formalism, and the `sys.enable` list.
- `L` — the Liouvillian used for time propagation; must be square and dimensionally consistent with `rho`.
- `rho` — the initial state vector for time propagation.
- `nstates` (optional) — if specified, the `nstates` most populated states and all mandatory unit coordinates are kept, irrespective of the tolerance parameter.

Output:
- `projector` — projector matrix into the reduced space (a column-subset of the identity, or the scalar `1` when elimination is skipped).

## References

- Source: [kernel/utilities/zte.m](https://github.com/IlyaKuprov/Spinach/blob/main/kernel/utilities/zte.m)
- Spinach Wiki: [zte.m](https://spindynamics.org/wiki/index.php?title=zte.m)
- I. Kuprov, zero track elimination method: [http://dx.doi.org/10.1016/j.jmr.2008.08.008](http://dx.doi.org/10.1016/j.jmr.2008.08.008)
