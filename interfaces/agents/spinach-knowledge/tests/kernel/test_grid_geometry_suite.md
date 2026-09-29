# tests/kernel/test_grid_geometry_suite.m

**Source:** [https://github.com/IlyaKuprov/Spinach/blob/main/tests/kernel/test_grid_geometry_suite.m](https://github.com/IlyaKuprov/Spinach/blob/main/tests/kernel/test_grid_geometry_suite.m)

## Purpose

Regression test suite for Spinach grid and spherical geometry helper functions. It verifies spherical arc and area formulae, Gauss-Legendre quadrature exactness, polar-grid structure, spherical quadrature weights, Voronoi solid angles, direct-product grid construction, SHREWD weights, and seeded repulsion-grid invariants.

## Behaviour

The function announces the test target with `fprintf`, creates a test result object via `new_test_result` under the name `kernel/grid_geometry_suite`, and then runs a sequence of assertions using `test_close` and `test_true`. Each assertion records a descriptive message. The checks performed are:

- **Spherical arc lengths** — `arclength(ex,ey)` must equal `pi/2` and `arclength(ex,-ex)` must equal `pi`, both to absolute and relative tolerances of `1e-15`, where `ex`, `ey`, `ez` are the Cartesian basis vectors on the unit sphere.
- **Spherical triangle area** — `sphtarea(ex,ey,ez,'unsigned')` must equal `pi/2` (the positive-octant triangle occupies one eighth of the sphere), and `sphtarea(ex,ez,ey,'signed')` must equal `-pi/2` (reversing vertex orientation flips the sign), both at tolerance `1e-15`.
- **Spherical triangle subdivision** — `sphtrsubd(ex,ey,ez)` returns three midpoints; the xy midpoint must equal `(ex+ey)/sqrt(2)` and all three subdivision points must have unit 2-norm, at tolerance `1e-15`.
- **Gauss-Legendre quadrature** — `gaussleg(-1,1,4)` must return all-positive weights, nodes sorted in ascending order, `sum(wg)` equal to 2, `sum(wg.*xg)` equal to 0, and `sum(wg.*xg.^2)` equal to `2/3`, with the integral checks at tolerance `1e-14`.
- **Polar grid** — `grid_polar(4,2)` must produce `3*(4-1)^2+2*(4-1)+1 = 28` points; the first radius must be 0 at tolerance `1e-15`; all radii must satisfy `r_pol>=0` and `max(r_pol)<=2`; the returned Laplacian must be of size `[npnts npnts]`; and `L_pol*ones(npnts,1)` must vanish at tolerance `1e-12` (constant function lies in the nullspace).
- **Fibonacci grid** — `grid_fibon('fib',3)` must return 7 points (the `fib` parameter `n` produces `2*n+1` points), zero alpha angles at tolerance `1e-15`, unit-norm Cartesian points at tolerance `1e-14`, all-positive Voronoi weights, and weights summing to 1 at tolerance `1e-12`.
- **Igloo grid** — `grid_igloo(4)` must return zero alpha angles at tolerance `1e-15`, unit vectors at tolerance `1e-14`, all-positive weights, and weights summing to 1 at tolerance `1e-12`.
- **Triangular grid** — `grid_trian('asg',2)` must return zero alpha angles at tolerance `1e-15` and unit vectors at tolerance `1e-14`; Voronoi weights are deliberately not computed for this grid because they are expensive.
- **Convex hull** — on a regular tetrahedron (vertices at `(±1,±1,±1)/sqrt(3)` with an even number of minus signs, converted to spherical angles via `acos` and `atan2`), `get_hull(theta,phi)` must return 4 hull facets and 12 directed edges (six undirected edges in both directions).
- **Voronoi solid angles** — `voronoisphere(xyz_tet)` must return all-positive solid angles summing to `4*pi` at tolerance `1e-12`, with four equal cells of `pi` steradians each by symmetry; `one_vcell_solidangle` on the first guided cell must match `sangles(1)`, and `vcell_solidangle(vertices,indices,xyz_tet)` must match the full vector of solid angles, both at tolerance `1e-12`.
- **Direct-product grid** — `grid_kron(angles1,weights1,angles2,weights2)` on two 2-point inputs must return 4 points, and the weights must equal `kron(weights1,weights2)` at tolerance `1e-15`.
- **SHREWD weights** — `shrewd(zeros(4,1),theta,phi,1,1e-8)` on the tetrahedral grid must return all-positive weights summing to 1 at tolerance `1e-12`.
- **Grid test** — `grid_test(0,0,0,1,0,'D_lmn')` must return a rank-zero residual of 0 at tolerance `1e-15` (a single unit-weight point integrates the constant Wigner function exactly).
- **Spatial grid point estimate** — `ngridpts(amps,durs,'1H',coh_order,sample_size)` with `amps=[0.01 -0.02]`, `durs=[1e-3 2e-3]`, `coh_order=2`, and `sample_size=0.01` must equal the reference value `ceil(abs(coh_order*spin('1H')*sum(abs(amps.*durs)))*sample_size/pi)` exactly (tolerance 0). The source comment states this is a change detector mirroring the implemented worst-case spiral formula, not an independent derivation.
- **Repulsion grid** — the global RNG state is saved with `rng`, reseeded with `rng(1,'twister')`, `repulsion(5,3,1)` is called, and the original RNG state is restored. The returned weights must equal `ones(5,1)/5` at tolerance `1e-15` (uniform weights), and the angles must describe unit vectors at tolerance `1e-14`.

## Inputs and outputs

```matlab
result=test_grid_geometry_suite()
```

- **Outputs**
  - `result` — regression test result object with explanatory messages, accumulated from the individual `test_close` and `test_true` assertions.
- **Inputs**
  - None.

## References

- [Spinach MATLAB source on GitHub — test_grid_geometry_suite.m](https://github.com/IlyaKuprov/Spinach/blob/main/tests/kernel/test_grid_geometry_suite.m)
- Functions exercised by the suite: `arclength`, `sphtarea`, `sphtrsubd`, `gaussleg`, `grid_polar`, `grid_fibon`, `grid_igloo`, `grid_trian`, `get_hull`, `voronoisphere`, `one_vcell_solidangle`, `vcell_solidangle`, `grid_kron`, `shrewd`, `grid_test`, `ngridpts`, `repulsion`, `spin`, `new_test_result`.
