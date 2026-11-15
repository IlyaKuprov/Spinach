# cpmg_echo

Syntax: `E=cpmg_echo(spin_system,parameters,H,R,K)`

Experiment owns validation, preparation and detection through cce_observation. times gives total elapsed sequence duration in seconds. Normalized signal defaults to true and uses the actual detected preparation reference. Explicit rho0/coil or preparation descriptors override the transition preparation. nloops is a positive integer; ideal transition-selective x pi pulses occur at (k-1/2)/nloops. Existing experiments/cpmg is untouched.
