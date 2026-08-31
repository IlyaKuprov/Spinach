# tests/kernel/test_grid_geometry_suite.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/tests/kernel/test_grid_geometry_suite.m`
- Signature: `result=test_grid_geometry_suite()`
- Total lines: 178

## Purpose

Tests grid and spherical geometry helpers. Syntax: result=test_grid_geometry_suite()

## Physical / mathematical content

- This file belongs to the `tests` part of Spinach. Its role should be read together with nearby files in the same directory, which usually share a common physical regime or infrastructure purpose.
- Orientation or trajectory averaging is performed numerically, so grid design, weights, and integration error control matter directly to accuracy and runtime.

## Numerical / algorithmic content

- Numerical integration over angles or geometry is part of the implementation, so point placement and weights are as important as the local Hamiltonian calculations.

## Code-derived implementation details

### Comment-guided execution stages

- Lines 17-18: Announce the test target; implemented by `fprintf('TESTING: Grid and spherical geometry helpers\n')`.
- Lines 20-23: State the grid target of the test; implemented by `result=new_test_result('kernel/grid_geometry_suite', 'Grid and spherical geometry helpers', 'grid helpers must return unit vectors, positive normalised weights, and exac…`.
- Lines 25-26: Define Cartesian basis vectors on the unit sphere; implemented by `ex=[1;0;0]`.
- Lines 30-32: Check spherical arc lengths between orthogonal and opposite points; implemented by `result=test_close(result,'arclength orthogonal axes',arclength(ex,ey),pi/2,1e-15,1e-15, 'orthogonal unit vectors subtend a quarter great circle')`.
- Lines 36-38: Check the area of the positive-octant spherical triangle; implemented by `result=test_close(result,'sphtarea octant unsigned',sphtarea(ex,ey,ez,'unsigned'),pi/2,1e-15,1e-15, 'the triangle bounded by three coordinate great circles occupies one…`.
- Lines 42-43: Check spherical midpoint subdivision; implemented by `[rxy,ryz,rzx]=sphtrsubd(ex,ey,ez)`.
- Lines 49-50: Check Gauss-Legendre quadrature on low-degree polynomials; implemented by `[xg,wg]=gaussleg(-1,1,4)`.
- Lines 62-63: Check polar grid point count, radial bounds, and Laplacian nullspace; implemented by `[~,r_pol,L_pol]=grid_polar(4,2)`.
- Lines 76-77: Check Fibonacci grid unit vectors and normalised Voronoi weights; implemented by `[alps,bets,gams,whts]=grid_fibon('fib',3)`.
- Lines 90-91: Check igloo grid weights on a small deterministic grid; implemented by `[alps,bets,gams,whts]=grid_igloo(4)`.
- Lines 102-103: Check a triangular grid without invoking the expensive Voronoi weights; implemented by `[alps,bets,gams]=grid_trian('asg',2)`.
- Lines 110-111: Check convex hull facets and directed edge list on a tetrahedron; implemented by `xyz_tet=[1 1 -1 -1;1 -1 1 -1;1 -1 -1 1]`.
- Lines 121-122: Check Voronoi solid angles on a regular tetrahedron; implemented by `[vertices,indices,~,sangles]=voronoisphere(xyz_tet)`.
- Lines 135-136: Check direct-product grid weights and cardinality; implemented by `angles1=[0 0 0;pi/2 0 0]`.
- Lines 146-147: Check SHREWD and grid-test invariants at rank zero and one; implemented by `weights=shrewd(zeros(4,1),theta,phi,1,1e-8)`.
- Lines 156-157: Check spatial-grid point estimate against the documented formula; implemented by `amps=[0.01 -0.02]`.
- Lines 165-166: Check seeded repulsion grid invariants without relying on a particular random cloud; implemented by `rng_state=rng`.

### Key state/data transformations

- Lines 21-23: computes `result` using `result=new_test_result('kernel/grid_geometry_suite', 'Grid and spherical geometry helpers', 'grid helpers must return unit vectors, positive normalised weights, and exac…`.
- Lines 26: computes `ex` using `ex=[1;0;0]`.
- Lines 27: computes `ey` using `ey=[0;1;0]`.
- Lines 28: computes `ez` using `ez=[0;0;1]`.
- Lines 43: computes `[rxy,ryz,rzx]` using `[rxy,ryz,rzx]=sphtrsubd(ex,ey,ez)`.
- Lines 50: computes `[xg,wg]` using `[xg,wg]=gaussleg(-1,1,4)`.
- Lines 63: computes `[~,r_pol,L_pol]` using `[~,r_pol,L_pol]=grid_polar(4,2)`.
- Lines 64: computes `npnts` using `npnts=3*(4-1)^2+2*(4-1)+1`.
- Lines 77: computes `[alps,bets,gams,whts]` using `[alps,bets,gams,whts]=grid_fibon('fib',3)`.
- Lines 78: computes `xyz` using `xyz=[sin(bets).*cos(gams) sin(bets).*sin(gams) cos(bets)]`.
- Lines 103: computes `[alps,bets,gams]` using `[alps,bets,gams]=grid_trian('asg',2)`.
- Lines 111: computes `xyz_tet` using `xyz_tet=[1 1 -1 -1;1 -1 1 -1;1 -1 -1 1]`.
- Lines 113: computes `theta` using `theta=acos(xyz_tet(3,:)).'`.
- Lines 114: computes `phi` using `phi=atan2(xyz_tet(2,:),xyz_tet(1,:)).'`.
- Lines 115: computes `[hull,edges]` using `[hull,edges]=get_hull(theta,phi)`.
- Lines 122: computes `[vertices,indices,~,sangles]` using `[vertices,indices,~,sangles]=voronoisphere(xyz_tet)`.
- Lines 129: computes `S_cell` using `S_cell=one_vcell_solidangle(vertices(:,indices{1}),xyz_tet(:,1))`.
- Lines 136: computes `angles1` using `angles1=[0 0 0;pi/2 0 0]`.

## Outputs

- result -regression test result with explanatory messages
- The test checks spherical arc and area formulae, Gauss-Legendre exactness,
- polar-grid structure, spherical quadrature weights, Voronoi solid angles,
- grid products, SHREWD weights, and seeded repulsion-grid invariants.

## Implementation structure

- Tests grid and spherical geometry helpers. Syntax:
- result=test_grid_geometry_suite()
- result -regression test result with explanatory messages
- The test checks spherical arc and area formulae, Gauss-Legendre exactness,
- polar-grid structure, spherical quadrature weights, Voronoi solid angles,
- grid products, SHREWD weights, and seeded repulsion-grid invariants.
- Announce the test target
- State the grid target of the test
- Define Cartesian basis vectors on the unit sphere
- Check spherical arc lengths between orthogonal and opposite points
- Check the area of the positive-octant spherical triangle
- Check spherical midpoint subdivision

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `new_test_result()`, `test_close()`, `arclength()`, `sphtarea()`, `sphtrsubd()`, `gaussleg()`, `test_true()`, `all()`, `diff()`, `grid_polar()`, `r_pol()`, `grid_fibon()`, `grid_igloo()`, `grid_trian()`, `acos()`, `xyz_tet()`.
