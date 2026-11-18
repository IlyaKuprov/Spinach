# test_clexpn_experiments

Syntax: `result=test_clexpn_experiments()`

Tests independent dense matrix propagation for free evolution, Hahn, CPMG and ESEEM, zero-time convention, bath pump boundaries, complex receiver phase and column-major reduced density.

Finite-drive regressions include the two-spin zero-drift z-drive acceptance probe, transverse preparation, Hilbert/Liouville equivalence, overlapping and partial intervals, detuning, fractional placement, noncommuting coupled drift and simultaneous GKSL relaxation. Unsupported ideal events and projected drives are rejected.

Preparation selectors verified for the example migration: `transition.ms=[-0.5 0.5]` for E4, `[-1 1]` for E3, `[-0.5 0.5;-0.5 0.5]` for two electrons, and `[-0.5 0.5;3.5 3.5]` for E/51V. `transition.spins` follows the observed-spin order. These select the equal two-manifold superpositions without manually assembled vectors.
Regression coverage includes native coherent error bounds and spectral residuals at zero tolerance, observed-transition parity against stock construction, lower-order numerical masks, and complete-bath validity when the requested error budget is deliberately smaller than the computed estimate.

Coherent projected Hahn spectral products are additionally compared with the general interval path for an arbitrary complex observed preparation, phased receiver and both normalization settings. Error auxiliaries must remain positive and finite.

The projected Hahn six-state electron/nitrogen bath checks zero-duration amplitude against the evaluated mixed-density trace with zero tolerance. This catches replacement by a mathematically equivalent literal one that changes the frozen floating-point result.

Named finite_drive coverage compares the declarative drive to the finite-control density bitwise and rejects invalid spins, axis, amplitude, detuning, duration, extra fields, early observation and competing events.

The finite-interval checks still supply raw sequence records through clexpn/free_evolution. That event-array interface is removed from the supported public contract, but its rejection is not enforced by these tests: they currently assert acceptance. Preserve the overlapping, fractional, partial-interval, coupled-drift, and dissipative physics checks when closing this remaining caller migration; a passing registered suite alone does not establish the no-legacy-field requirement.
