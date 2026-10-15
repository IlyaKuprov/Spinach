# Remaining deterministic kernel utility tests

`tests/kernel/test_dynamic_remaining_core_suite.m` returns a regression result
for deterministic utility contracts that are not covered by more specialised
suites. Failed assertions are retained in the result rather than silently
accepted.

The physical checks cover adiabatic elimination and its induced relaxation,
interaction-frame transformations, dipolar geometry and units, isotope
rescaling, and magnetic pumping. The Clebsch–Gordan check uses the shared
`clebsch_gordan` routine and verifies the positive spin-half triplet coefficient
for opposite projections in the documented coupled-state-first convention.

Numerical checks cover finite-difference Hessians, Gaussian and Lorentzian
line shapes (including degenerate geometries and narrow widths), relaxation
kite selection, frozen-state columns, the Sorensen bound, and trajectory
stitching. Random-walk checks verify finite Euler angles and an identity
initial orientation. Metadata, reporting, cell packaging, and polyadic
diagnostics are exercised alongside those numerical contracts. Synthetic
compiled fixtures use per-substance descriptors and offsets.

The source test retains the exact fixtures and tolerances; this entry records
the guarantees rather than duplicating that executable specification.
