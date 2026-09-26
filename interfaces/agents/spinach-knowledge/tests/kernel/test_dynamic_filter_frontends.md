# tests/kernel/test_dynamic_filter_frontends.m

- Signature: `result=test_dynamic_filter_frontends()`

## Purpose

Tests dynamic state-filter front-end kernels on compact spin systems. Syntax: result=test_dynamic_filter_frontends()

## Physical / mathematical content

## Numerical / algorithmic content

- The file is built around the standard Spinach workflow: create the spin system, choose a basis or context, assemble operators/superoperators, then propagate or analyse the resulting dynamics.

## Outputs

- result -regression test result with explanatory messages
- The test exercises coherence(), correlation(), decouple(), homospoil(),
- and spinlock() on small spherical-tensor Liouville-space systems with
- analytically known surviving state components.

## Implementation structure

- Tests dynamic state-filter front-end kernels on compact spin systems. Syntax:
- result=test_dynamic_filter_frontends()
- result -regression test result with explanatory messages
- The test exercises coherence(), correlation(), decouple(), homospoil(),
- and spinlock() on small spherical-tensor Liouville-space systems with
- analytically known surviving state components.
- Announce the test target
- State the dynamic filter target of the test
- Check analytical coherence and correlation filters
- Check analytical decoupling of a coupled spin system
- Check homospoil zero-quantum and longitudinal filters
- Check analytical spin-lock projection along both transverse axes
