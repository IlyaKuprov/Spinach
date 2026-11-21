# cpmg_echo

Syntax: `[E,aux]=cpmg_echo(spin_system,parameters,H,R,K)`

Experiment owns validation, preparation, and detection through cce_observation. times gives total elapsed sequence duration in seconds. Normalised signal defaults to true and uses the actual detected preparation reference. Explicit rho0/coil or preparation descriptors override the transition preparation. nloops is a positive integer; ideal transition-selective x pi pulses occur at (k-1/2)/nloops. Existing experiments/cpmg is untouched.

Caller-supplied parameters.sequence is rejected, including empty arrays. Use finite_drive, refocusing_train, or the appropriate named echo experiment; internal events are constructed from native physical declarations, never read from that removed field.

Detection uses the fixed physical receiver: after an odd observed-spin manifold interchange, the propagated bath-coherence block is adjointed, including its initial complex preparation coefficient. The receiver itself and the initial normalisation reference are not conjugated. Even nloops retains the original coherence ordering; odd nloops requires the adjoint. Zero elapsed time applies no pulse.

The declared signal-normalisation default is resolved before the local input checks and physical preparation.

`aux.abs_error_bound` sums the floating-point state propagation estimates over the free intervals. Ideal inversions preserve the Frobenius norm. Each receiver scales the sum by its Frobenius norm and, when normalised, by the inverse magnitude of its initial reference. Zero elapsed time contributes no propagation error; the estimate is not a certified error bound.
