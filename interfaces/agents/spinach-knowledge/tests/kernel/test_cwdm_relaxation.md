# tests/kernel/test_cwdm_relaxation.m

Compares two distinct heteronuclear substances with CSA and dipolar couplings
against independently constructed single-substance spin systems. Complete
relaxation matrices agree within 1e-12 relative Frobenius norm for Redfield
(kite and secular retention, asynchronous and serial integration), Lindblad,
SRFK, and combined T1/T2 plus Redfield. The substances have different isotopes,
interaction tensors, Redfield correlation times, and phenomenological rates.
The same fixture checks the Hamiltonian direct sum and exactly zero unit rows
and columns in every unthermalised relaxation block.

The test records that `weak` is not a supported relaxation theory. Nottingham
with two electrons in one substance and a nucleus in another, or with the
electrons split between substances, must raise the documented substance-local
error. An ordinary two-electron-plus-nucleus substance must still produce
nonzero relaxation with a conserved unit coordinate. This tests an explicit
unsupported-domain boundary, not a Nottingham direct-sum numerical equality.
Compiled descriptors with two electron pairs also exercise the rejection,
with contiguous and interleaved global electron ordering. These are built
without a relaxation theory, then given the supported single-substance
Nottingham relaxation settings to test the consumer independently of
`create`'s existing two-electron restriction.

A two-substance IME T1/T2 fixture checks both `steady` methods against
independently constructed single-substance equilibria, with empty and nonzero
initial guesses. Every unit coordinate stays at its specified
chemical concentration. Later-block trace-row and initial-normalisation
violations are rejected. Both solvers must also reject an identity-propagator
block paired with a thermalised block, with each substance tested in turn.

Pumping targets in the second substance and both substances are checked
against their own unit columns. Unequal unit populations distinguish local
pumping from an accidental source through the first substance. Identity
components in a later block must be rejected.

NGCE rejects the segmented two-substance system with both zero and nonzero
regularisation using `Spinach:ngce:segmentedSubstances`. A supported
single-substance zero stochastic trajectory retains zero uncertainty and
regularises only the non-unit directions.

Solid-effect coverage asserts the named segmented steady-state boundary with an electron and an independent nucleus. A supported single-substance zero-drive calculation is compared with thermal-equilibrium detection.

A two-substance Hilbert fixture converted by `sim2liouv` checks the named segmented Zeeman steady boundary for both methods, with empty and supplied initial guesses.

Both DNP scan functions are called with complete two-substance inputs and must raise their named segmented-substance boundary rather than attempting a singular single-trace solve.

Steady-state comparisons and guesses carry the specified concentrations, including the spin-free pool. Pumping targets use unweighted coil_state; their source columns act on the instantaneous populations.

A symmetric relaxation fixture with cross-substance active-state entries and exactly zero unit action must raise `Spinach:thermalize:crossSubstanceRelaxation` in IME. The corresponding supported block-diagonal operator must annihilate the supplied polarised unweighted target within `10*eps*norm(R_therm,'fro')*norm(rho_eq)` in the whole-vector two-norm.

Pumping reference shapes use the explicit four-argument `coil_state` API.
