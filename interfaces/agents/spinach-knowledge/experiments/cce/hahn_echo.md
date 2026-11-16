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

For generalized Zeeman-Hilbert dynamics with zero R and K and no bath pump, the experiment diagonalizes the Hermitian retained-cluster Hamiltonian once, composes the two free intervals with the physical transition pi rotation, and retains observed densities in aux.density. aux also reports absolute numerical bounds, eigensystem and orthogonality residuals, the Hamiltonian assembly bound and backend. The registered evaluator retains this evidence in raw cluster records. Bounds use the same spectral and arithmetic operation count as the legacy generalized Hahn evaluator; density output is column-vectorized and scalar/multiple detectors retain their normalization divisors.
