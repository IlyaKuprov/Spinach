# examples/microfluidics/plain_diff.m

Source: [examples/microfluidics/plain_diff.m](https://github.com/IlyaKuprov/Spinach/blob/main/examples/microfluidics/plain_diff.m)

## System and transport model

This example transports longitudinal proton magnetisation on a cropped COMSOL chip mesh. The mesh and velocity files are `chip_mesh.txt` and `chip_velo.txt`; after import, both mesh-velocity components are set to zero, leaving diffusion and distal-pipe drainage. The example is marked as having no dynamics in the spin subspace and has no chemical reaction network.

The mesh import uses crop ranges `[286.8 287.5]` and `[576.0 579.0]` and the source explicit inactive-cell mask. The propagation arrays use `2659` cells. The source does not state coordinate units. A spatially uniform Hamiltonian and relaxation operator are supplied across cells; the signal is represented by the Lz state.

## Initial condition and parameters

The initial Lz magnetisation is `0.5` in cells `1240` and `1246`, zero elsewhere. The coil detects Lz uniformly in all cells. The time step is `50` seconds and the trajectory has `500` points. Diffusion is `1e-7 m^2/s`. Drainage is assigned as `-0.01` from cell `2600` through the final cell; the source does not annotate the drainage-rate unit. It is applied through `K_op={speye([4 4])}` and the per-cell `K_ph` profile.

## Calculation and observable

`meshflow(spin_system,@simple_flow,parameters)` propagates the spatial trajectory. `fpl2phan` reshapes the Lz coil signal to a `[2659 parameters.npoints]` cell-by-time array. Each time point is rendered as a mesh concentration map labelled in arbitrary units; the script animates these maps rather than saving a spectrum or numeric result file.

## Scope

With velocity explicitly zeroed, the simulated transport is diffusion plus distal-pipe drainage, without advection or chemical conversion.
