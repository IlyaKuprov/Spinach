# cce_eseem

Syntax: `E=cce_eseem(spin_system,parameters,H,R,K)`

Experiment owns validation, preparation and detection through cce_observation. times gives total elapsed sequence duration in seconds. Normalized signal defaults to true and uses the actual detected preparation reference. Explicit rho0/coil or preparation descriptors override the transition preparation. Two-pulse ESEEM with the equivalent transition preparation and a midpoint refocusing pulse. The name avoids collision with the existing eseem experiment.
