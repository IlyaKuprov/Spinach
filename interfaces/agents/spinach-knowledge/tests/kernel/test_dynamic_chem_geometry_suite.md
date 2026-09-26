# tests/kernel/test_dynamic_chem_geometry_suite.m

- Signature: `result=test_dynamic_chem_geometry_suite()`

## Purpose

Tests deterministic chemistry and geometry utility helpers. Syntax: result=test_dynamic_chem_geometry_suite()

## Physical / mathematical content

## Numerical / algorithmic content

## Outputs

- result -regression test result with explanatory messages
- The test checks lattice construction, geometric measurements, coupling
- extraction, chemical shifts, nearest-neighbour lookup, and tensor helpers.

## Implementation structure

- Tests deterministic chemistry and geometry utility helpers. Syntax:
- result=test_dynamic_chem_geometry_suite()
- result -regression test result with explanatory messages
- The test checks lattice construction, geometric measurements, coupling
- extraction, chemical shifts, nearest-neighbour lookup, and tensor helpers.
- Announce the test target
- State the chemistry and geometry target of the test
- Build a small spin-system descriptor for metadata helpers
- Check simple cubic lattice construction and periodic vectors
- Check a signed right-angle dihedral from four Cartesian points
- Check Cartesian point-cloud binning on a regular grid
- Check nearest-neighbour lookup by Cartesian coordinate distance
