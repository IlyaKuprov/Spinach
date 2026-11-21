# cce_eseem

Syntax: `[E,aux]=cce_eseem(spin_system,parameters,H,R,K)`

Experiment owns validation, preparation, and detection through cce_observation. times gives total elapsed sequence duration in seconds. Normalised signal defaults to true and uses the actual detected preparation reference. Explicit rho0/coil or preparation descriptors override the transition preparation. Two-pulse ESEEM with the equivalent transition preparation and a midpoint refocusing pulse. The name avoids collision with the existing eseem experiment.

Caller-supplied parameters.sequence is rejected, including empty arrays. Use finite_drive, refocusing_train, or the appropriate named echo experiment; internal events are constructed from native physical declarations, never read from that removed field.

Detection uses the fixed physical receiver: after an odd observed-spin manifold interchange, the propagated bath-coherence block is adjointed, including its initial complex preparation coefficient. The receiver itself and the initial normalisation reference are not conjugated. The single midpoint inversion requires this adjoint at positive elapsed time. Zero elapsed time applies no pulse.

The declared signal-normalisation default is resolved before the local input checks and physical preparation.

`aux.abs_error_bound` sums the floating-point state propagation estimates over the free intervals. Ideal inversions preserve the Frobenius norm. Each receiver scales the sum by its Frobenius norm and, when normalised, by the inverse magnitude of its initial reference. Zero elapsed time contributes no propagation error; the estimate is not a certified error bound.
