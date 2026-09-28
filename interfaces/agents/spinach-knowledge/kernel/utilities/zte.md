# kernel/utilities/zte.m

- Signature: `projector=zte(spin_system,L,rho,nstates)`

## Purpose

Zero track elimination inspects the first few steps of a trajectory and removes states that remain below the drop tolerance. If `nstates` is specified, it instead keeps the states with the greatest trajectory weight.

## Inputs and output

- `L`: Liouvillian used for time propagation.
- `rho`: initial state vector used for time propagation.
- `nstates` (optional): number of states to keep, irrespective of the tolerance-based selection. It must be a positive integer no greater than the state-space dimension.
- `projector`: projector into the reduced space. Use `L_reduced=P'*L*P` and `rho_reduced=P'*rho`, with `P=projector`.

`L` and `rho` must be numeric; `L` must be square, and `rho` must be a single column state vector with matching dimension. Zero track elimination is available only for the `zeeman-liouv` and `sphten-liouv` formalisms.

## Algorithm and notes

The function chooses a time step of `1/cheap_norm(L)` (or a unit time step if this is infinite), then propagates up to `spin_system.tols.zte_nsteps` trajectory points using the Krylov-based `step` function. It stops early when the number of states whose maximum sampled amplitude exceeds `spin_system.tols.zte_tol` stops changing. Without `nstates`, it drops states whose maximum sampled amplitude is below that tolerance; with `nstates`, it retains the specified number of states with the largest maximum sampled amplitudes.

The default tolerance may be changed by setting `sys.tols.zte_tol` before calling `create.m`. The basis is left unchanged if zero track elimination is disabled, the state vector has too few zeros under `spin_system.tols.zte_maxden`, or its 1-norm is below `spin_system.tols.zte_tol`. If tiny interactions or nearly equivalent spins are present, disable zero track elimination by adding `'zte'` to `sys.disable`.

Further information is available in IK's JMR paper: http://dx.doi.org/10.1016/j.jmr.2008.08.008. See also <https://spindynamics.org/wiki/index.php?title=zte.m>.