# cce_acquire — native acquisition service

Syntax: `answer=cce_acquire(spin_system,parameters,H,R,K)`

Shared preparation, propagation, and detection for free_evolution, finite_drive, and refocusing_train. Native physical declarations are compiled by cce_pulse_events; caller-supplied sequence fields are rejected, including empty arrays. times gives total elapsed sequence duration in seconds. Normalized signal defaults to true and uses the actual detected preparation reference. Explicit rho0/coil or preparation descriptors override the transition preparation. Direct density output is time-by-column-vectorised observed density; the context uses its native density-reduction contract.

GKSL relaxation follows `L=H+1i*R+1i*K`. Projected propagation acts on an alpha/beta coherence block, not a physical density. Generalized density output is column-vectorised and reduced elementwise, including the empty reference; positivity is not guaranteed at truncated order.

Finite_drive and refocusing_train own physical control declarations. This
service compiles internal events from those validated native declarations. Hilbert propagation
uses cce_segments for retained interval Hamiltonians, endpoint rotations, and
exterior Cartesian-field updates. Its auxiliary density and segment arrays are
available to direct experiment callers; the context retains numerical evidence,
not an old public result wrapper. Overlapping finite amplitudes add, and partial
pulses end at the requested sample. The existing pure-finite Liouville path
retains drift and GKSL relaxation. Projected driven propagation is refused.


`preparation=transverse` means the positive trace-one density `(identity+2*Lx)/2` for exactly one observed spin one-half; it is not a deviation density. Other populations use the existing explicit state descriptor or expert `rho0`. Detection still needs a transition, a descriptor or an expert `coil`; density output supplies the identity detector automatically.
## Coherent numerical diagnostics

The optional second output in `[answer,aux]=cce_acquire(spin_system,parameters,H,R,K)` carries `abs_error_bound`, `eigen_residual`, `orthogonality_residual`, `hamiltonian_bound` and `backend` for conventional coherent free evolution. One eigensystem is reused for all requested times; each state is evaluated from the initial density, never recursively. Residuals use the original Hamiltonian, preserving real versus complex arithmetic. Finite-drive and projected branches retain their existing propagation.

## Generalized Liouville density transport

For density output without finite events, the zeeman-liouv path uses the stock step exponential action (or stock propagator when explicitly selected), with prop_chop=0. Its auxiliary absolute error estimate is 64*eps*dimension*(1+norm(generator,1)*time) per density element. This preserves the established master-equation numerical path and its diagnostics. Coherent Hilbert, projected, controlled, and signal-only paths retain their existing contracts.

Controlled Hilbert evolution uses segmentwise matrix exponentials in the original multiplication order, including the initial zero-duration interval. Its existing arithmetic estimate is 16*d^3*eps times the sum of 1+abs(dt)*norm(H,2) over segments; density entries inherit that estimate. The named finite_drive experiment compiles its declarative rectangular drive through this path.

For an uncontrolled, nondissipative projected transition with default preparation and detection, the experiment diagonalizes each conditional bath Hamiltonian once. It propagates the two branches in their original product order and returns `trace(rho_bath)` exactly at zero elapsed time. The auxiliary eigenvalue residuals, orthogonality residuals and forward bounds use the retained bath dimension. Explicit preparation/detection overrides and controlled or dissipative propagation retain their general paths.

Nonprojected propagation is reduced to the observed-spin density before detection, preserving the partial-trace and detector multiplication order. The Liouville-space action and its floating-point error estimate apply to signal and observable outputs as well as density output. Observable bounds use the observed detector norm and the actual normalization reference.

Coherent non-projected Liouville free evolution supplies a structurally diagonal generator in sparse storage to stock `step`. Explicit master-equation inputs retain dense generator arithmetic, including diagonal dephasing. Projected master evolution uses the same stock action and 64*eps*n*(1+norm(L,1)*time) arithmetic estimate as generalized master evolution, with the actual detection scale. Normalized projected preparation is divided by its reference before propagation; unnormalized preparation retains its amplitude. The explicit `propagator=stock` choice remains supported.

Explicit multi-observable finite drives compile their detector cells before propagation. Error-norm factors follow output-column order for row and column receiver cells.

Projected selected-transition FIDs use the same spectral manifold trace even with an explicit observed preparation. Explicit zeeman-liouv selects the master-equation action also when the jump list is empty; normalized projected propagation starts from the bath density and detects its trace. Coherent generalized master generators use sparse stock action, while active dissipative generators remain dense. This refines the earlier diagonal-only storage description: sparsity follows the absence of bath jumps, not isdiag alone.
