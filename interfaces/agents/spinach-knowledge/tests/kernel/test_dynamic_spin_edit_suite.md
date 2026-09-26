# tests/kernel/test_dynamic_spin_edit_suite.m

- Signature: `result=test_dynamic_spin_edit_suite()`

## Purpose

Tests deterministic spin-system editing support utilities. Syntax: result=test_dynamic_spin_edit_suite()

## Physical / mathematical content

- The relaxation model is Redfield-type perturbation theory: fluctuating interactions enter through correlation functions or spectral densities and generate a linear relaxation superoperator.

## Numerical / algorithmic content

## Outputs

- result -regression test result with explanatory messages
- The test checks spin removal, dilute-isotope subsystem generation,
- assumption overrides, and merging of Spinach input structures.

## Implementation structure

- Tests deterministic spin-system editing support utilities. Syntax:
- result=test_dynamic_spin_edit_suite()
- result -regression test result with explanatory messages
- The test checks spin removal, dilute-isotope subsystem generation,
- assumption overrides, and merging of Spinach input structures.
- Announce the test target
- State the spin editing target of the test
- Build a small but structurally complete spin-system descriptor
- Check numeric spin removal updates identities, labels, coordinates, and parts
- Check chemical part renumbering across multiple subsystems
- Check destruction of stale basis, connectivity, symmetry, and assumption data
- Check logical spin removal follows the same path
