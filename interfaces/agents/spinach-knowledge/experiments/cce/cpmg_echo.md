# cpmg_echo

Syntax: `E=cpmg_echo(spin_system,parameters,H,R,K)`

Experiment owns validation, preparation and detection through cce_observation. times gives total elapsed sequence duration in seconds. Normalized signal defaults to true and uses the actual detected preparation reference. Explicit rho0/coil or preparation descriptors override the transition preparation. nloops is a positive integer; ideal transition-selective x pi pulses occur at (k-1/2)/nloops. Existing experiments/cpmg is untouched.

Caller-supplied parameters.sequence is rejected, including empty arrays. Use finite_drive, refocusing_train, or the appropriate named echo experiment; internal events are constructed from native physical declarations, never read from that removed field.

Projected evolution restores the physical detected coherence by taking the bath-coherence adjoint after an odd inversion count. Even counts remain in the original ordering, and zero elapsed time applies no pulse. Complex preparation and receiver phase therefore retain the same convention as full-space evolution.
