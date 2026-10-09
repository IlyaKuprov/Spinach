# tests/kernel/test_cwdm_thermalisation.m

T5 and T6 evolve two independent proton pairs at 14.1 T and 298 K under their Hamiltonian and IME relaxation. Unit-coordinate errors at five times through 20 s and the final equilibrium infinity-norm error must be below 1e-10. A third spin-free block checks inert population storage. The generator must be exactly independent of the initial concentrations, including zero populations.

The public `thermalize` call is checked against `relaxation` using unit-concentration thermal targets. Passing preweighted targets must raise `Spinach:thermalize:targetConcentration`; nonunit columns must remain unchanged.

T7 explicitly collapses the independent identity coordinates into a common identity and combines the unit-concentration thermal sources. Hamiltonian and nonunit dissipative intertwining residuals must be exactly zero, the unit-column residual must be nonzero, and the longitudinal recovery ratio must be within 1e-5 of two for this near-identical proton-pair fixture. This algebraic representation test does not call a second installation or establish a stock-kernel comparison; that measurement is performed separately.

The single-substance Zeeman-Liouville test also rejects target traces zero, 0.3, and two, while a unit-trace target makes the concentration-weighted equilibrium stationary under isotropic damping.

Detection uses the explicit `exact` method of `coil_state`.
