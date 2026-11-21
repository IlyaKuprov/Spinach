# cce_acquire — native acquisition service

[MATLAB source](../../../../../kernel/cce/cce_acquire.m)

Syntax: `answer=cce_acquire(spin_system,parameters,H,R,K)`

Shared preparation, propagation, and detection for free_evolution, finite_drive, and refocusing_train. Native physical declarations are compiled by cce_pulse_events; caller-supplied sequence fields are rejected, including empty arrays. times gives total elapsed sequence duration in seconds. Normalised signal defaults to true and uses the actual detected preparation reference. Explicit rho0/coil or preparation descriptors override the transition preparation. Direct density output is time-by-column-vectorised observed density; the context uses its native density-reduction contract.

GKSL relaxation follows `L=H+1i*R+1i*K`. Projected propagation acts on an alpha/beta coherence block, not a physical density. Generalised density output is column-vectorised and reduced elementwise, including the empty reference; positivity is not guaranteed at truncated order.

Finite_drive and refocusing_train own physical control declarations. This
service compiles internal events from those validated native declarations. Hilbert propagation
uses cce_segments for retained interval Hamiltonians, endpoint rotations, and
exterior Cartesian-field updates. Its auxiliary density and segment arrays are
available to direct experiment callers; the context retains numerical evidence,
not an old public result wrapper. Overlapping finite amplitudes add, and partial
pulses end at the requested sample. The existing pure-finite Liouville path
retains drift and GKSL relaxation. Projected driven propagation is refused.

`preparation=transverse` means the positive trace-one density `(identity+2*Lx)/2` for exactly one observed spin one-half; it is not a deviation density. Other populations use the existing explicit state descriptor or expert `rho0`. Detection still needs a transition, a descriptor, or an expert `coil`; density output supplies the identity detector automatically.

## Coherent numerical diagnostics

The optional second output in `[answer,aux]=cce_acquire(spin_system,parameters,H,R,K)` carries `abs_error_bound`, `eigen_residual`, `orthogonality_residual`, `hamiltonian_bound`, and `backend` for conventional coherent free evolution. One eigensystem is reused for all requested times; each state is evaluated from the initial density, never recursively. Residuals use the original Hamiltonian, preserving real versus complex arithmetic. Finite-drive and projected branches retain their existing propagation.

## Generalised Liouville density transport

For density output without finite events, the zeeman-liouv path uses the stock step exponential action (or stock propagator when explicitly selected), with prop_chop=0. Its auxiliary absolute error estimate is 64*eps*dimension*(1+norm(generator,1)*time) per density element. This preserves the established master-equation numerical path and its diagnostics. Coherent Hilbert, projected, controlled, and signal-only paths retain their existing contracts.

Controlled Hilbert evolution uses segmentwise matrix exponentials in the original multiplication order, including the initial zero-duration interval. Its existing arithmetic estimate is 16*d^3*eps times the sum of 1+abs(dt)*norm(H,2) over segments; density entries inherit that estimate. The named finite_drive experiment compiles its declarative rectangular drive through this path.

For an uncontrolled, nondissipative projected transition with selected-transition detection, the experiment diagonalises each conditional bath Hamiltonian once, including when the observed preparation is explicit. It propagates the two branches in their original product order; normalised output returns `trace(rho_bath)` exactly at zero elapsed time, while unnormalised output retains the detected preparation reference. The auxiliary eigenvalue residuals, orthogonality residuals, and forward bounds use the retained bath dimension. Explicit detector overrides and controlled or dissipative propagation use their general paths.

For a coherent projected FID with an explicit transition-compatible coil, the propagated state is the vectorised bath-coherence block. Its Frobenius-norm error estimate is multiplied by the Frobenius norm of the trace detector (including its complex coefficient) and divided by the magnitude of the applied normalisation divisor. This propagation-only estimate is zero at zero elapsed time; it is not a rigorous bound on the total signal error.

Nonprojected propagation is reduced to the observed-spin density before detection, preserving the partial-trace and detector multiplication order. The Liouville-space action and its floating-point error estimate apply to signal and observable outputs as well as density output. Observable bounds use the observed detector norm and the actual normalisation reference.

Coherent non-projected Liouville free evolution without bath jumps supplies its generator in sparse storage to stock `step`. Nonempty bath-jump declarations retain dense generator arithmetic, including diagonal dephasing. Projected master evolution uses the same stock action and 64*eps*n*(1+norm(L,1)*time) arithmetic estimate as generalised master evolution, with the actual detection scale. Normalised projected preparation is divided by its reference before propagation; unnormalised preparation retains its amplitude. The explicit `propagator=stock` choice remains supported.

Explicit multi-observable finite drives compile their detector cells before propagation. Error-norm factors follow output-column order for row and column receiver cells.

Projected selected-transition FIDs use the same spectral manifold trace even with an explicit observed preparation. Explicit zeeman-liouv selects the master-equation action also when the jump list is empty; normalised projected propagation starts from the bath density and detects its trace. Coherent generalised master generators use sparse stock action when bath-jump declarations are absent or empty; nonempty declarations retain dense arithmetic. Storage selection does not depend on `isdiag`.

Finite Liouville drives retain dense segmentwise exponential propagation and accumulate the floating-point estimate `64*eps*n*sum(1+abs(dt)*norm(L_segment,1))` over nonzero intervals, including free drift before and after each pulse. Here n is the retained Liouville dimension; overlapping drives enter the same segment generator. Density entries use this estimate directly; signals and observables apply the observed detector norm and the actual normalisation reference. Zero elapsed time requires no segment propagation. The estimate is marked `floating_point_estimate` and reaches raw records, full-bath reduction, and numerical-budget diagnostics; it is not a rigorous truncation bound.
