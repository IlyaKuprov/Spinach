# hahn_echo

Syntax: `E=hahn_echo(spin_system,parameters,H,R,K)`

Experiment owns validation, preparation and detection through cce_observation. times gives total elapsed sequence duration in seconds. Normalized signal defaults to true and uses the actual detected preparation reference. Explicit rho0/coil or preparation descriptors override the transition preparation. One ideal transition-selective x pi pulse at half the total time.

GKSL relaxation follows `L=H+1i*R+1i*K`. Projected propagation acts on an alpha/beta coherence block, not a physical density. A branch-swapping pulse exchanges both Hamiltonians and the projected dissipator. Detection and normalization remain experiment-owned.

For an odd number of observed-spin branch swaps, the projected trace labels the opposite coherence from a fixed generalized receiver. Preserve the explicitly selected receiver convention; complex traces may be conjugates even when magnitudes agree.

Hahn echo also supports generalized `output='density'`, with column-vectorised observed density rows and the same elementwise reduction/full-bath override as free evolution. Projected coherence blocks do not expose a density.

The optional second output carries spectral residuals, orthogonality residuals and absolute numerical error estimates for coherent projected Hahn echoes. The two invariant-branch spectra are resolved once per retained cluster. Branch propagator products preserve the physical ideal-pulse phase and the zero-time identity; no legacy driver or event-array adapter is called.

Projected zero-duration Hahn propagation evaluates the retained bath density trace, preserving floating-point identity for mixed baths whose dimension is not a power of two. The six-state electron/nitrogen regression checks this exactly.

The polarized two-proton example supplies sampled Zeeman-basis populations through parameters.cce.bath_state and averages after expansion. The returned signal remains separate from estimator standard error and deterministic summation bounds in diagnostics.reduction.

## Generalized unitary Hahn propagation

By default, for generalized Zeeman-Hilbert dynamics with zero R and K and no bath pump, the experiment diagonalizes the Hermitian retained-cluster Hamiltonian once, composes the two free intervals with the physical transition pi rotation, and retains observed densities in aux.density. aux also reports absolute numerical bounds, eigensystem and orthogonality residuals, the Hamiltonian assembly bound and backend. The registered evaluator retains this evidence in raw cluster records. Bounds use the same spectral and arithmetic operation count as the legacy generalized Hahn evaluator; density output is column-vectorized and scalar/multiple detectors retain their normalization divisors.

For default selected-transition detection in the generalized coherent path, the response is evaluated as `alpha'*rho_observed*beta`, with the receiver phase and the detected preparation reference applied once. Its detector norm is one for the normalized transition states. Explicit coils and detection descriptors retain the Hilbert–Schmidt trace evaluation and their own operator norms.

For generalized Hilbert evolution, parameters.cce.numerics.propagator='expm' selects direct segment exponentials with the same physical midpoint pi pulse. The initial zero-duration interval is retained in the arithmetic bound 16*d^3*eps*sum(1+abs(dt)*norm(H,2)). The default spectral route and its residual-derived bounds are unchanged; direct-segment output does not invent spectral residuals.

The direct-segment selection rejects projected, dissipative, non-Hilbert, and pump-bearing requests rather than ignoring that numerical preference.

Native ideal master-equation echoes use stock step or the explicitly selected stock propagator on each interval, retaining sparse coherent generalized generators and dense dissipative generators. They return observed densities separately from normalized signals and report the original dimension/norm arithmetic estimate. A scalar refocusing descriptor replaces the midpoint observed pi rotation in generalized evolution: spins identifies the complete observed system or one physical site, operator is its Hermitian integrated generator, and fraction places it on the total-time axis. This preserves selective off-midpoint bath rotations without exposing an event-array interface; correlated observed/bath operator pulses and projected selective rotations are rejected. Capability version 4 invalidates older Hahn propagation checkpoints.

The refocusing descriptor also accepts an absolute time in seconds instead of fraction. At zero elapsed time no pulse acts; nonzero samples include a pulse exactly at their endpoint and omit later pulses. Projected pulses are accepted only when their physical matrix preserves or interchanges the selected manifolds without leakage, extending the earlier generalized-only selective contract. aux.pulse_time and aux.pulse_applied retain the actual sample-specific placement and inclusion decisions. Hahn capability version 5 supersedes version 4.

Independent full-density exponentials also check absolute refocusing before, at, and after its timestamp for both eig and expm generalized propagation. The direct interval path applies its rotation only when that sample includes the pulse.
