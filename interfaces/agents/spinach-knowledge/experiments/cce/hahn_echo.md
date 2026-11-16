# hahn_echo

Syntax: `E=hahn_echo(spin_system,parameters,H,R,K)`

Experiment owns validation, preparation and detection through cce_observation. times gives total elapsed sequence duration in seconds. Normalized signal defaults to true and uses the actual detected preparation reference. Explicit rho0/coil or preparation descriptors override the transition preparation. One ideal transition-selective x pi pulse at half the total time.

GKSL relaxation follows `L=H+1i*R+1i*K`. Projected propagation acts on an alpha/beta coherence block, not a physical density. A branch-swapping pulse exchanges both Hamiltonians and the projected dissipator. Detection and normalization remain experiment-owned.

For an odd number of observed-spin branch swaps, the projected trace labels the opposite coherence from a fixed generalized receiver. Preserve the explicitly selected receiver convention; complex traces may be conjugates even when magnitudes agree.

Hahn echo also supports generalized `output='density'`, with column-vectorised observed density rows and the same elementwise reduction/full-bath override as free evolution. Projected coherence blocks do not expose a density.
