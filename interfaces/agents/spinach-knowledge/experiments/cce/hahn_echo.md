# hahn_echo

[MATLAB source](../../../../../experiments/cce/hahn_echo.m)

Syntax: `E=hahn_echo(spin_system,parameters,H,R,K)`

Experiment owns validation, preparation, and detection through cce_observation. times gives total elapsed sequence duration in seconds. Normalised signal defaults to true and uses the actual detected preparation reference. Explicit rho0/coil or preparation descriptors override the transition preparation. One ideal transition-selective x pi pulse at half the total time.

GKSL relaxation follows `L=H+1i*R+1i*K`. Projected propagation acts on an alpha/beta coherence block, not a physical density. A branch-swapping pulse exchanges both Hamiltonians and the projected dissipator. Detection and normalisation remain experiment-owned.

Detection uses the fixed physical receiver: after an odd observed-spin manifold interchange, the propagated bath-coherence block is adjointed, including its initial complex preparation coefficient. The receiver itself and the initial normalisation reference are not conjugated. Zero elapsed time applies no inversion; manifold-preserving pulses retain the original coherence ordering.

Hahn echo also supports generalised `output='density'`, with column-vectorised observed density rows and the same elementwise reduction/full-bath override as free evolution. Projected coherence blocks do not expose a density.

The optional second output carries spectral residuals, orthogonality residuals, and absolute numerical error estimates for coherent projected Hahn echoes. The two invariant-branch spectra are resolved once per retained cluster. Branch propagator products preserve the physical ideal-pulse phase and the zero-time identity; no legacy driver or event-array adapter is called.

Projected zero-duration Hahn propagation evaluates the retained bath density trace, preserving floating-point identity for mixed baths whose dimension is not a power of two. The six-state electron/nitrogen regression checks this exactly.

The polarised two-proton example supplies sampled Zeeman-basis populations through parameters.bath_probs with parameters.bath_samples and averages after expansion. The returned signal remains separate from estimator standard error and deterministic summation bounds in diagnostics.reduction.

## Generalised unitary Hahn propagation

By default, for generalised Zeeman-Hilbert dynamics with zero R and K and no bath pump, the experiment diagonalises the Hermitian retained-cluster Hamiltonian once, composes the two free intervals with the physical transition pi rotation, and retains observed densities in aux.density. aux also reports absolute numerical bounds, eigensystem and orthogonality residuals, the Hamiltonian assembly bound, and backend. The registered evaluator retains this evidence in raw cluster records. Bounds use the same spectral and arithmetic operation count as the legacy generalised Hahn evaluator; density output is column-vectorised and scalar/multiple detectors retain their normalisation divisors.

For default selected-transition detection in the generalised coherent path, the response is evaluated as `alpha'*rho_observed*beta`, with the receiver phase and the detected preparation reference applied once. Its detector norm is one for the normalised transition states. Explicit coils and detection descriptors retain the Hilbert–Schmidt trace evaluation and their own operator norms.

For generalised Hilbert evolution, parameters.cce.propagator='expm' selects direct segment exponentials with the same physical midpoint pi pulse. The initial zero-duration interval is retained in the arithmetic bound 16*d^3*eps*sum(1+abs(dt)*norm(H,2)). The default spectral route and its residual-derived bounds are unchanged; direct-segment output does not invent spectral residuals.

The direct-segment selection rejects projected, dissipative, non-Hilbert, and pump-bearing requests rather than ignoring that numerical preference.

Native ideal master-equation echoes use stock step or the explicitly selected stock propagator on each interval, retaining sparse coherent generalised generators and dense dissipative generators. They return observed densities separately from normalised signals and report the original dimension/norm arithmetic estimate. A scalar refocusing descriptor replaces the midpoint observed pi rotation in generalised evolution: spins identifies the complete observed system or one physical site, operator is its Hermitian integrated generator, and fraction places it on the total-time axis. This preserves selective off-midpoint bath rotations without exposing an event-array interface; correlated observed/bath operator pulses and projected selective rotations are rejected.

The refocusing descriptor also accepts an absolute time in seconds instead of fraction. At zero elapsed time no pulse acts; nonzero samples include a pulse exactly at their endpoint and omit later pulses. Projected pulses are accepted only when their physical matrix preserves or interchanges the selected manifolds without leakage, extending the earlier generalised-only selective contract. aux.pulse_time and aux.pulse_applied retain the actual sample-specific placement and inclusion decisions.

Independent full-density exponentials also check absolute refocusing before, at, and after its timestamp for both eig and expm generalised propagation. The direct interval path applies its rotation only when that sample includes the pulse.

Refocusing may alternatively specify axis and angle for a Cartesian physical-spin rotation instead of an integrated operator. A one-site exterior Cartesian rotation updates its product density and the subsequent static mean field in generalised Hilbert dynamics; direct segment exponentials preserve interval timing and arithmetic bounds. Exterior explicit-operator pulses remain rejected, so arbitrary correlated exterior operations are not approximated by frozen fields.

Caller-supplied parameters.sequence is rejected, including empty arrays. Use finite_drive, refocusing_train, or the appropriate named echo experiment; internal events are constructed from native physical declarations, never read from that removed field.

Checkpoint capability versions distinguish the fixed physical receiver; regenerate records produced under the preceding spin-manifold-overlap convention.

Projected coherent, generalised Hilbert, and master-equation propagation are explicit branches in the experiment body. The physical pulse remains fixed across direct segment propagation; receiver conventions and arithmetic bounds are unchanged.

The interval fallback accumulates free-evolution and two-sided density-rotation roundoff. An exact projected branch exchange costs no arithmetic; zero elapsed time skips the pulse. Each observable uses its retained receiver Frobenius norm and actual normalisation divisor. For a reduced density, each partial-trace component instead has sensitivity `sqrt(bath_dimension)`. These are propagation estimates, not preparation or readout error bounds.
