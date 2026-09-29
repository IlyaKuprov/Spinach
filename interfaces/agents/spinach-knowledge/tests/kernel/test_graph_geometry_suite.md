# tests/kernel/test_graph_geometry_suite.m

**Source:** [https://github.com/IlyaKuprov/Spinach/blob/main/tests/kernel/test_graph_geometry_suite.m](https://github.com/IlyaKuprov/Spinach/blob/main/tests/kernel/test_graph_geometry_suite.m)

## Purpose

Regression test for graph, geometry, lattice, and coordinate utilities in Spinach. It verifies small graph decompositions, lattice generation, dihedral angles, coordinate density binning, nearest-spin lookup, substance and label lookup helpers, and coordinate-derived dipolar tensor identities.

## Behaviour

The function announces the test target with `fprintf`, initialises a test result via `new_test_result` with the identifier `kernel/graph_geometry_suite` and the description "Graph, geometry, and coordinate utilities", then runs a sequence of assertions:

- **Cubic lattice:** calls `cubic_lattice('13C',1.5,2)` and checks that the isotope count is 8 (a two-period cubic lattice contains 2^3 isotope entries), that the sorted coordinates enumerate all corners of the requested cubic grid (compared against `1.5*[0 0 0;0 0 1;0 1 0;0 1 1;1 0 0;1 0 1;1 1 0;1 1 1]` with tolerances `1e-15`), and that the first periodic boundary vector equals `[3 0 0]` (spacing times number of periods along x).
- **Dihedral angle:** computes `dihedral([1 0 0],[0 0 0],[0 1 0],[0 1 1])` and checks that the magnitude of the angle is 90 degrees (tolerance `1e-13`), since orthogonal adjacent planes give a ninety-degree dihedral magnitude.
- **Depth-first partitioning:** applies `dfpt` with size 2 to the three-node path graph `logical([0 1 0;1 0 1;0 1 0])` and checks the result equals `logical([1 1 0;0 1 1])` — the two edges of the path.
- **Strongly connected components:** applies `scomponents` to two disconnected two-cycles (`logical([0 1 0 0;1 0 0 0;0 0 0 1;0 0 1 0])`) and checks that nodes 1 and 2 share one component, nodes 3 and 4 share another, and the two components differ.
- **Point-density binning:** calls `xyz2pd` on three points (`[0.25 0.25 0.25;0.75 0.75 0.75;1.25 0.25 0.25]`) over the unit cube with a 2×2×2 grid and checks the counts match a reference where cell (1,1,1) and cell (2,2,2) each contain one point and the point outside the grid is discarded (tolerance `1e-15`).
- **Nearest-spin lookup:** in the three-spin geometry fixture, `nearest_spin(spin_system,1)` returns index 3 at distance 0.5 — the Euclidean coordinate separation to the closest coordinate.
- **Label and substance lookup:** checks that `idxof(sys,'ca')` returns 2 for a system with labels `{'ha','ca','hb'}` and isotopes `{'1H','13C','1H'}`, and that `which_subst(spin_system,[1 3])` returns 1 because spins one and three belong to the first chemical substance.
- **Dipolar tensor identities:** calls `xyz2dd([0 0 0],[0 0 1],'1H','1H')` and checks that the Euler angles alpha, beta, gamma are all zero for a z-axis dipolar vector (tolerance `1e-13`), that the dipolar matrix equals the coupling times `diag([1 1 -2])` (the axial traceless dipolar tensor, tolerances `1e-10`/`1e-12`), and that doubling the internuclear distance to `[0 0 2]` divides the dipolar coupling by eight (tolerance `1e-8`/`1e-12`).
- **Hyperfine tensor:** calls `xyz2hfc([0 0 0],[0 0 1],'1H')` and checks that the point-dipole hyperfine tensor is symmetric (tolerance `1e-12`) and traceless (tolerance `1e-12`).

The three-spin geometry fixture contains spins (`'1H'`, `'13C'`, `'1H'`), chemical parts `{[1 3],2}`, and coordinates `{[0 0 0],[2 0 0],[0.5 0 0]}`.

## Inputs and outputs

```matlab
result=test_graph_geometry_suite()
```

- **Output:** `result` — regression test result object with explanatory messages, accumulated through the `test_true` and `test_close` assertions.
- **Input:** none.

## References

- [Spinach MATLAB source on GitHub](https://github.com/IlyaKuprov/Spinach/blob/main/tests/kernel/test_graph_geometry_suite.m)
