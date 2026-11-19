# test_clexpn_experiments

[MATLAB source](../../../../../tests/kernel/test_clexpn_experiments.m)

Syntax: `result=test_clexpn_experiments()`

Tests independent dense matrix propagation for free evolution, Hahn, CPMG, and ESEEM, zero-time convention, bath pump boundaries, complex receiver phase, and column-major reduced density.

Finite-drive regressions include the two-spin zero-drift z-drive acceptance probe, transverse preparation, Hilbert/Liouville equivalence, overlapping and partial intervals, detuning, fractional placement, noncommuting coupled drift, and simultaneous GKSL relaxation. Raw finite, ideal, and empty event arrays and projected finite drives are rejected.

Preparation selectors verified for the example migration: `transition.ms=[-0.5 0.5]` for E4, `[-1 1]` for E3, `[-0.5 0.5;-0.5 0.5]` for two electrons, and `[-0.5 0.5;3.5 3.5]` for E/51V. `transition.spins` follows the observed-spin order. These select the equal two-manifold superpositions without manually assembled vectors.
Regression coverage includes native coherent error bounds and spectral residuals at zero tolerance, observed-transition parity against stock construction, lower-order numerical masks, and complete-bath validity when the requested error budget is deliberately smaller than the computed estimate.

Coherent projected Hahn spectral products are additionally compared with the general interval path for an arbitrary complex observed preparation, phased receiver, and both normalisation settings. Error auxiliaries must remain positive and finite.

The projected Hahn six-state electron/nitrogen bath checks zero-duration amplitude against the evaluated mixed-density trace with zero tolerance. This catches replacement by a mathematically equivalent literal one that changes the frozen floating-point result.

Named finite_drive coverage compares the declarative drive to the finite-control density bitwise and rejects invalid spins, axis, amplitude, detuning, duration, extra fields, and competing events, and checks partial observation against independent coupled propagation.

Every finite-interval physical check uses the native finite_drive experiment, including independent-axis overlaps, fractional starts, partial observations, coupled drift, and simultaneous GKSL relaxation. Explicit rejection assertions cover removed raw finite, ideal, and empty sequence inputs.

Complex transverse preparation and a phased receiver compare projected CPMG/ESEEM with full-space evolution for odd and even inversion counts in Hilbert and Liouville representations, including zero elapsed time.

Fixed-receiver Hahn regressions compare complete complex amplitudes with independent full-density or generalised propagation. Coverage includes a bath and an empty bath, complex observed preparation, Hilbert and Liouville representations, spectral and interval paths, and normalised and unnormalised detection; zero time retains the initial reference.

Nonzero bath relaxation checks projected master echoes with manifold-interchanging and preserving rotations. Bath pumping constructs its rotation in Hilbert space even when the state is propagated in Liouville space.

CPMG and ESEEM full-bath regressions compare each receiver’s accumulated estimate to an independent unitary norm calculation. They distinguish zero elapsed time from positive-time numerical budget failure in both Zeeman representations and in projected and generalised propagation.

Non-spectral Hahn and bath-pumped Hahn tests check strict full-bath budgets, coincident and projected pulse arithmetic, complex scaled receivers, and the partial-trace scale for density output. An independent amplitude-damping calculation checks Frobenius gain and the two interval norms for projected and generalised bath-pumped echoes.
