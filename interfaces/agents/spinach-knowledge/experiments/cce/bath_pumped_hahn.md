# bath_pumped_hahn

Syntax: `E=bath_pumped_hahn(spin_system,parameters,H,R,K)`

Experiment owns validation, preparation and detection through cce_observation. times gives total elapsed sequence duration in seconds. Normalized signal defaults to true and uses the actual detected preparation reference. Explicit rho0/coil or preparation descriptors override the transition preparation. Hahn refocusing plus an ideal bath x pi pulse at pump.fraction. pump requires center, halfwidth, units (Hz or rad/s), fraction and allow_empty. Optional sites/species restrict the physical bath. Spectral boundaries are inclusive. Exterior conditioned-spin pumping is rejected.
