# tests/kernel/test_graph_geometry_suite.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/tests/kernel/test_graph_geometry_suite.m`
- Signature: `result=test_graph_geometry_suite()`
- Total lines: 120

## Purpose

Tests graph, geometry, lattice, and coordinate utilities. Syntax: result=test_graph_geometry_suite()

## Physical / mathematical content

- This file belongs to the `tests` part of Spinach. Its role should be read together with nearby files in the same directory, which usually share a common physical regime or infrastructure purpose.
- The spin physics includes through-space magnetic dipole-dipole coupling, a rank-2 anisotropic interaction with strong orientation dependence and characteristic secular/non-secular structure.

## Numerical / algorithmic content

- The file also defines local helper function(s): `local_geometry_system()`. This usually means the public entry point is supported by tightly coupled validation or helper logic kept private to the file.

## Outputs

- result -regression test result with explanatory messages
- The test checks small graph decompositions, lattice generation,
- dihedral angles, coordinate density binning, nearest-spin lookup,
- substance and label lookup helpers, and coordinate-derived dipolar
- tensor identities.

## Implementation structure

- Tests graph, geometry, lattice, and coordinate utilities. Syntax:
- result=test_graph_geometry_suite()
- result -regression test result with explanatory messages
- The test checks small graph decompositions, lattice generation,
- dihedral angles, coordinate density binning, nearest-spin lookup,
- substance and label lookup helpers, and coordinate-derived dipolar
- tensor identities.
- Announce the test target
- State the utility target of the test
- Check a two-period cubic lattice coordinate and periodic-cell layout
- Check a right-handed coordinate set with a ninety-degree dihedral
- Check depth-first partitioning on a three-node path graph

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `new_test_result()`, `cubic_lattice()`, `sortrows()`, `cell2mat()`, `test_true()`, `test_close()`, `dihedral()`, `logical()`, `dfpt()`, `isequal()`, `scomponents()`, `component_idx()`, `xyz2pd()`, `density_ref()`, `local_geometry_system()`, `nearest_spin()`.
