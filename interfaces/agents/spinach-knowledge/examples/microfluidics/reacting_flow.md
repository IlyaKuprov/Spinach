# examples/microfluidics/reacting_flow.m

## Model and transport

This example evolves concentrations on a cropped, inactivated COMSOL mesh with imported velocity data; it has no spin dynamics. Two competing second-order reaction channels are coupled to the shared flow-diffusion generator, with a fifth inert solvent component. The source sets `k1=2.0` toward exo and `k2=1.0` toward endo (bimolecular rate units `L/(mol*s)` for mol/L concentrations and seconds) and diffusion to `1e-7` (no unit is stated for that value). The mesh crop is x=`[286.8, 287.5]` and y=`[576.0, 579.0]`; the code also inactivates listed mesh indices. It seeds component 1 at cell 1240 with 0.50 and component 2 at cell 1246 with 0.25; these initial values have no unit specified, and the plotted fields are labelled a.u.

## Time stepping and observable

The script takes 280 steps of 20 seconds. At each step it evaluates the concentration-dependent reaction matrix in every cell, combines those matrices with the flow-diffusion generator, and advances the flattened state with `step`. Four mesh concentration fields are animated in a 2-by-2 plot: cyclopentadiene, acrylonitrile, exo-NBCN, and endo-NBCN. The source header describes runtime as seconds; this is not a measured runtime. No separate boundary-condition rule is stated beyond using the imported, cropped and inactivated mesh.

## Source

https://github.com/IlyaKuprov/Spinach/blob/main/examples/microfluidics/reacting_flow.m
