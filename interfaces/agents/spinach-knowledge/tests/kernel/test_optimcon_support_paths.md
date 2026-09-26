# tests/kernel/test_optimcon_support_paths.m

- Signature: `result=test_optimcon_support_paths()`

## Purpose

Tests small optimal-control support paths. Syntax: result=test_optimcon_support_paths()

## Physical / mathematical content

## Numerical / algorithmic content

- Time propagation is explicit. In Spinach this usually means repeated application of matrix exponentials or propagator factorizations to density operators or state vectors in Hilbert/Liouville/Fokker-Planck space.
- Finite-difference discretisation appears in the implementation, so numerical accuracy depends on stencil order, boundary handling, and the balance between resolution and conditioning.

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
