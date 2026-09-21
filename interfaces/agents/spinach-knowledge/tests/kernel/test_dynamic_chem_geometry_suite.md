# tests/kernel/test_dynamic_chem_geometry_suite.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/tests/kernel/test_dynamic_chem_geometry_suite.m`
- Signature: `result=test_dynamic_chem_geometry_suite()`
- Total lines: 130

## Purpose

Tests deterministic chemistry and geometry utility helpers. Syntax: result=test_dynamic_chem_geometry_suite()

## Physical / mathematical content

- This file belongs to the `tests` part of Spinach. Its role should be read together with nearby files in the same directory, which usually share a common physical regime or infrastructure purpose.

## Numerical / algorithmic content

- The file also defines local helper function(s): `local_geometry_spin_system()`. This usually means the public entry point is supported by tightly coupled validation or helper logic kept private to the file.

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

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `new_test_result()`, `local_geometry_spin_system()`, `cubic_lattice()`, `test_true()`, `all()`, `strcmp()`, `test_close()`, `isequal()`, `dihedral()`, `xyz2pd()`, `density_ref()`, `nearest_spin()`, `which_subst()`, `get_coupling()`, `chemshifts()`, `cs_ppm()`.
