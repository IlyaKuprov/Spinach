# examples/microfluidics/plain_flow.m

Source: [examples/microfluidics/plain_flow.m](https://github.com/IlyaKuprov/Spinach/blob/main/examples/microfluidics/plain_flow.m)

## System and transport model

This example transports longitudinal proton magnetisation on a COMSOL-imported microfluidic chip mesh with its velocity field. The source describes the case as having no dynamics in the spin subspace. It adds diffusion and a distal-pipe drainage sink, but no chemical reaction network. The crop ranges are `[286.8 287.5]` and `[576.0 579.0]`; the import also applies the explicit inactive-cell list in the source. The coordinate units are not stated.

The model uses a common Hamiltonian and relaxation operator across the `2659` cells, via spatial profiles of ones. It represents the signal as Lz, not as a spatially resolved set of chemical species.

## Initial condition and parameters

The initial Lz profile is `2` in cells `140:160`, zero elsewhere; detection uses Lz in all cells. The time step is `50` seconds and `npoints=200`. Diffusion is `1e-7 m^2/s`. The drainage profile is zero except for `-0.01` at cells `2600:end`, applied with `K_op={speye([4 4])}` and the per-cell `K_ph` profile. The source gives no unit for this drainage value.

## Calculation and observable

`meshflow(spin_system,@simple_flow,parameters)` computes the trajectory using the imported flow field and the specified diffusion and drainage. `fpl2phan` reshapes the detected Lz signal to `[2659 parameters.npoints]`. The script animates `real(traj(:,n))` as a mesh concentration map labelled in arbitrary units; it does not save a separate data file.

## Scope

The example transports one injected Lz profile through a prescribed COMSOL velocity field; chemical conversion is outside this model.
