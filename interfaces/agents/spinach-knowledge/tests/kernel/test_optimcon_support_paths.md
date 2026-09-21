# tests/kernel/test_optimcon_support_paths.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/tests/kernel/test_optimcon_support_paths.m`
- Signature: `result=test_optimcon_support_paths()`
- Total lines: 352

## Purpose

Tests small optimal-control support paths. Syntax: result=test_optimcon_support_paths()

## Physical / mathematical content

- This file belongs to the `tests` part of Spinach. Its role should be read together with nearby files in the same directory, which usually share a common physical regime or infrastructure purpose.

## Numerical / algorithmic content

- Time propagation is explicit. In Spinach this usually means repeated application of matrix exponentials or propagator factorizations to density operators or state vectors in Hilbert/Liouville/Fokker-Planck space.
- Finite-difference discretisation appears in the implementation, so numerical accuracy depends on stencil order, boundary handling, and the balance between resolution and conditioning.
- The file also defines local helper function(s): `local_spin_system()`, `local_data()`, `local_penalty_term()`, `local_finite_grad()`, `local_finite_hess()`, `local_penalty_grad()`, `local_objective()`, `local_total_objective()`. This usually means the public entry point is supported by tightly coupled validation or helper logic kept private to the file.

## Outputs

- result -regression test result with explanatory messages
- The test covers penalty functions, trapezium-product derivatives,
- objective-function collection, and small line-search helper paths.

## Implementation structure

- Tests small optimal-control support paths. Syntax:
- result=test_optimcon_support_paths()
- result -regression test result with explanatory messages
- The test covers penalty functions, trapezium-product derivatives,
- objective-function collection, and small line-search helper paths.
- Announce the test target
- State the support-path target of the test
- Make a minimal quiet Spinach object for low-level helper calls
- Check the no-penalty path
- Check the norm-square penalty against its closed form
- Check the spillout penalty against explicit clipping residuals
- Check the derivative norm-square gradient by finite differences

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `new_test_result()`, `penalty()`, `trapdiff()`, `objeval()`, `local_spin_system()`, `test_close()`, `waveform()`, `spill_hi()`, `spill_lo()`, `spill_mask()`, `local_finite_grad()`, `local_penalty_term()`, `test_true()`, `grad_dns()`, `isequal()`, `amp_waveform()`.
