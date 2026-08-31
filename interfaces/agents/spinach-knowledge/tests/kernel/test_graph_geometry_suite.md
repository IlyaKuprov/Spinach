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

## Code-derived implementation details

### Comment-guided execution stages

- Lines 18-19: Announce the test target; implemented by `fprintf('TESTING: Graph, geometry, and coordinate utilities\n')`.
- Lines 21-24: State the utility target of the test; implemented by `result=new_test_result('kernel/graph_geometry_suite', 'Graph, geometry, and coordinate utilities', 'Small graph and geometry helpers must reproduce exact combinatorial a…`.
- Lines 26-27: Check a two-period cubic lattice coordinate and periodic-cell layout; implemented by `[sys,inter]=cubic_lattice('13C',1.5,2)`.
- Lines 39-40: Check a right-handed coordinate set with a ninety-degree dihedral; implemented by `phi=dihedral([1 0 0],[0 0 0],[0 1 0],[0 1 1])`.
- Lines 45-46: Check depth-first partitioning on a three-node path graph; implemented by `path_graph=logical([0 1 0;1 0 1;0 1 0])`.
- Lines 52-53: Check strongly connected components on two disconnected two-cycles; implemented by `strong_graph=logical([0 1 0 0;1 0 0 0;0 0 0 1;0 0 1 0])`.
- Lines 61-62: Check three-dimensional point-density binning with one discarded point; implemented by `points=[0.25 0.25 0.25;0.75 0.75 0.75;1.25 0.25 0.25]`.
- Lines 71-72: Check nearest-spin lookup from Cartesian coordinates; implemented by `spin_system=local_geometry_system()`.
- Lines 79-80: Check label and chemical-substance lookup helpers; implemented by `sys.labels={'ha','ca','hb'}`.
- Lines 88-89: Check coordinate-derived dipolar coupling tensor identities on the z axis; implemented by `[dip_coupling,alpha,beta,gamma,dip_mat]=xyz2dd([0 0 0],[0 0 1],'1H','1H')`.
- Lines 101-102: Check coordinate-derived hyperfine tensor symmetry and tracelessness; implemented by `hfc_one=xyz2hfc([0 0 0],[0 0 1],'1H')`.

### Key state/data transformations

- Lines 22-24: computes `result` using `result=new_test_result('kernel/graph_geometry_suite', 'Graph, geometry, and coordinate utilities', 'Small graph and geometry helpers must reproduce exact combinatorial a…`.
- Lines 27: computes `[sys,inter]` using `[sys,inter]=cubic_lattice('13C',1.5,2)`.
- Lines 28: computes `coord_obs` using `coord_obs=sortrows(cell2mat(inter.coordinates))`.
- Lines 29: computes `coord_ref` using `coord_ref=sortrows(1.5*[0 0 0;0 0 1;0 1 0;0 1 1;1 0 0;1 0 1;1 1 0;1 1 1])`.
- Lines 40: computes `phi` using `phi=dihedral([1 0 0],[0 0 0],[0 1 0],[0 1 1])`.
- Lines 46: computes `path_graph` using `path_graph=logical([0 1 0;1 0 1;0 1 0])`.
- Lines 47: computes `subgraphs` using `subgraphs=dfpt(path_graph,2)`.
- Lines 48: computes `subgraph_ref` using `subgraph_ref=logical([1 1 0;0 1 1])`.
- Lines 53: computes `strong_graph` using `strong_graph=logical([0 1 0 0;1 0 0 0;0 0 0 1;0 0 1 0])`.
- Lines 54: computes `component_idx` using `component_idx=scomponents(strong_graph)`.
- Lines 62: computes `points` using `points=[0.25 0.25 0.25;0.75 0.75 0.75;1.25 0.25 0.25]`.
- Lines 63: computes `density_obs` using `density_obs=xyz2pd(points,[0 1],[0 1],[0 1],2,2,2)`.
- Lines 64: computes `density_ref` using `density_ref=zeros(2,2,2)`.
- Lines 65: computes `density_ref(1,1,1)` using `density_ref(1,1,1)=1`.
- Lines 66: computes `density_ref(2,2,2)` using `density_ref(2,2,2)=1`.
- Lines 72: computes `spin_system` using `spin_system=local_geometry_system()`.
- Lines 73: computes `[near_idx,near_dist]` using `[near_idx,near_dist]=nearest_spin(spin_system,1)`.
- Lines 80: computes `sys.labels` using `sys.labels={'ha','ca','hb'}`.

### Local helper functions

- Line 111: `local_geometry_system()` — `function spin_system=local_geometry_system()`.
  - Representative operation: `spin_system.comp.nspins=3`.
  - Representative operation: `spin_system.comp.isotopes={'1H','13C','1H'}`.

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
