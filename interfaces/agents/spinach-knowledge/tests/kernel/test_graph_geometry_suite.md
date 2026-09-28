# tests/kernel/test_graph_geometry_suite.m

- Signature: `result=test_graph_geometry_suite()`

## Purpose

Checks graph, lattice, geometry, coordinate, and chemical-label/substance helpers against small reference cases.

## Physical / mathematical content

From Cartesian coordinates, the suite checks a z-axis point-dipole coupling tensor and the symmetry and tracelessness of the point-dipole hyperfine tensor.

## Numerical / algorithmic content

The suite checks lattice coordinates and periodic-cell layout, dihedral-angle conventions, graph partitions, coordinate-density bins, nearest-spin lookup, and chemical label/substance lookup. Tensor checks compare with axial, symmetric, and traceless reference identities.

## Outputs

`result` is the regression-test result with explanatory messages, covering small graph decompositions, lattice generation, dihedral angles, coordinate-density binning, nearest-spin lookup, substance and label lookup helpers, and coordinate-derived dipolar tensor identities.

## Implementation structure

The tests progress through exact geometry and graph examples, coordinate-based lookup helpers, and tensor identities. Graph coverage includes depth-first partitioning on a three-node path and strongly connected components on two disconnected two-cycles.
