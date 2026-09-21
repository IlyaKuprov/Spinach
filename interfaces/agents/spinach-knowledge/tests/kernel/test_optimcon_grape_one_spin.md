# tests/kernel/test_optimcon_grape_one_spin.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/tests/kernel/test_optimcon_grape_one_spin.m`
- Signature: `result=test_optimcon_grape_one_spin()`
- Total lines: 147

## Purpose

Tests one-spin optimal-control setup and Hilbert-space GRAPE. Syntax: result=test_optimcon_grape_one_spin()

## Physical / mathematical content

- This file belongs to the `tests` part of Spinach. Its role should be read together with nearby files in the same directory, which usually share a common physical regime or infrastructure purpose.
- The numerical method is quasi-Newton optimisation: curvature information is approximated from successive step and gradient differences instead of forming exact second derivatives every iteration.
- The numerical method is limited-memory quasi-Newton optimisation, which keeps only a short curvature history and is therefore suitable for waveform vectors too large for dense Hessians.
- The optimisation logic is Newton or Newton-like: search directions use first- and second-order local curvature information, usually with regularisation or line-search safeguards.
- The control theory content is GRAPE: fidelity derivatives are propagated through a piecewise-constant pulse sequence so that waveform samples can be improved by gradient-based optimisation.

## Numerical / algorithmic content

- Time propagation is explicit. In Spinach this usually means repeated application of matrix exponentials or propagator factorizations to density operators or state vectors in Hilbert/Liouville/Fokker-Planck space.
- Finite-difference discretisation appears in the implementation, so numerical accuracy depends on stencil order, boundary handling, and the balance between resolution and conditioning.
- The file also defines local helper function(s): `local_spin_system()`. This usually means the public entry point is supported by tightly coupled validation or helper logic kept private to the file.

## Outputs

- result -regression test result with explanatory messages
- The test checks that optimcon() accepts a minimal one-spin Hilbert-space
- control problem, and that grape_hilb() fidelity and gradient agree with
- independent matrix exponentiation and finite differences.

## Implementation structure

- Tests one-spin optimal-control setup and Hilbert-space GRAPE. Syntax:
- result=test_optimcon_grape_one_spin()
- result -regression test result with explanatory messages
- The test checks that optimcon() accepts a minimal one-spin Hilbert-space
- control problem, and that grape_hilb() fidelity and gradient agree with
- independent matrix exponentiation and finite differences.
- Announce the test target
- State the optimal-control target of the test
- Ensure that a parallel pool is available for the ensemble loop
- Build a minimal Hilbert-space Spinach object for setup and GRAPE
- Define one-spin operators and a two-step timing grid
- Configure a minimal optimal-control problem

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `new_test_result()`, `optimcon()`, `grape_hilb()`, `gcp()`, `parpool()`, `local_spin_system()`, `pauli()`, `test_true()`, `test_close()`, `rmfield()`, `strcmp()`, `contains()`, `waveform()`, `pulse_dt()`, `hdot()`, `wf_plus()`.
