# tests/kernel/test_graph_geometry_suite.m

- Signature: `result=test_graph_geometry_suite()`

## Purpose

Tests graph, geometry, lattice, and coordinate utilities. Syntax: result=test_graph_geometry_suite()

## Physical / mathematical content

- The spin physics includes through-space magnetic dipole-dipole coupling, a rank-2 anisotropic interaction with strong orientation dependence and characteristic secular/non-secular structure.

## Numerical / algorithmic content

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
