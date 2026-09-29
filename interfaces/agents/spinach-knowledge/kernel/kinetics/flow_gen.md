# kernel/kinetics/flow_gen.m

- Signature: `F=flow_gen(spin_system,parameters)`
- Direct MATLAB source: [`kernel/kinetics/flow_gen.m`](https://github.com/IlyaKuprov/Spinach/blob/main/kernel/kinetics/flow_gen.m)
- Existing Wiki: [`flow_gen.m`](https://spindynamics.org/wiki/index.php?title=flow_gen.m)

## Purpose

Builds a sparse spatial-motion generator on the Voronoi mesh in `spin_system.mesh`, combining mesh advection and optional diffusion.

## Inputs and domain

`spin_system.mesh` must contain `idx` and `vor`; the implementation also uses active-cell and triangle indices, Voronoi vertices/cells/weights, mesh coordinates `x,y`, and velocity components `u,v`. `parameters.diff` is a finite, non-negative, real scalar diffusion coefficient in m^2/s; it defaults to zero. The source's `grumble` guard checks the mesh/index/Voronoi fields and diffusion value.

## Transfer calculation and units

For cell `k` and a neighbour `m` sharing exactly two Voronoi vertices, let `A_k=mesh.vor.weights(k)`, `b_km` be the shared-edge length, `r_km=(x_m-x_k,y_m-y_k)`, and `vbar_km=(v_k+v_m)/2`. The code computes the signed advection contribution

`q_km = -(b_km/(2 A_k ||r_km||)) * dot(vbar_km,r_km)`.

If `q_km>0`, it adds `+q_km` at matrix entry (k,m); otherwise it adds `-q_km` at (m,k). Diffusion is added at entry (m,k) as

`d_km = (b_km/(A_k ||r_km||)) * parameters.diff`.

With lengths in metres, velocity in m/s, and diffusion in m^2/s, both contributions have units s^-1. This is a spatial generator, not a frequency-domain line shape; the routine does not convert Hz to angular frequency or vice versa.

## Assembly and output

The directed contributions are assembled into an `ncells×ncells` sparse matrix `B`. The code subtracts each column sum on the diagonal, then applies Voronoi-weight scaling:

`F = diag(1./weights) * (B - diag(sum(B,1))) * diag(weights)`.

The returned `F` has one row and column per Voronoi cell. No normalisation of the mesh weights or state vector is performed here; the displayed diagonal balance and left/right weight scaling are the generator's actual construction.

## Construction guards

Only pairs sharing exactly two Voronoi vertices enter the local transfer calculation. Missing mesh, indexing, or Voronoi structures and an invalid diffusion value raise source-defined errors. Mesh geometry and field consistency beyond these explicit checks are not asserted here.
