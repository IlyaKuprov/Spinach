# free_evolution

Syntax: `answer=free_evolution(spin_system,parameters,H,R,K)`

Experiment owns validation, preparation and detection through cce_observation. times gives total elapsed sequence duration in seconds. Normalized signal defaults to true and uses the actual detected preparation reference. Explicit rho0/coil or preparation descriptors override the transition preparation. No controls. Direct density output is time-by-column-vectorised observed density; density context reduction is not yet available.

GKSL relaxation follows `L=H+1i*R+1i*K`. Projected propagation acts on an alpha/beta coherence block, not a physical density. Generalized density output is column-vectorised and reduced elementwise, including the empty reference; positivity is not guaranteed at truncated order.

Finite controls are owned by this experiment. `parameters.sequence` is an optional array of finite events: `type=finite`, physical `target` IDs, `axis` x/y/z or a unit three-vector, `amplitude` (rad/s), optional `detuning` (rad/s, default zero), `placement` absolute/fraction, `at` (seconds/fraction), and positive `duration` (seconds). Overlapping drives add; intervals are clipped at the requested elapsed time, and drift and GKSL relaxation remain active. Both generalized Zeeman representations are supported. Projected drives, ideal events and targets outside the retained subsystem are rejected rather than silently ignored. Named echo experiments retain ownership of their ideal controls.

`preparation=transverse` means the positive trace-one density `(identity+2*Lx)/2` for exactly one observed spin one-half; it is not a deviation density. Other populations use the existing explicit state descriptor or expert `rho0`. Detection still needs a transition, a descriptor or an expert `coil`; density output supplies the identity detector automatically.
## Coherent numerical diagnostics

The optional second output in `[answer,aux]=free_evolution(spin_system,parameters,H,R,K)` carries `abs_error_bound`, `eigen_residual`, `orthogonality_residual`, `hamiltonian_bound` and `backend` for conventional coherent free evolution. One eigensystem is reused for all requested times; each state is evaluated from the initial density, never recursively. Residuals use the original Hamiltonian, preserving real versus complex arithmetic. Finite-drive and projected branches retain their existing propagation.

## Generalized Liouville density transport

For density output without finite events, the zeeman-liouv path uses the stock step exponential action (or stock propagator when explicitly selected), with prop_chop=0. Its auxiliary absolute error estimate is 64*eps*dimension*(1+norm(generator,1)*time) per density element. This preserves the established master-equation numerical path and its diagnostics. Coherent Hilbert, projected, controlled, and signal-only paths retain their existing contracts.
