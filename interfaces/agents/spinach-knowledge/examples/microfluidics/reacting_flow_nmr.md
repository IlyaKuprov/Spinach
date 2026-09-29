# examples/microfluidics/reacting_flow_nmr.m

## Coupled model

A spatial Diels–Alder reaction model on imported COMSOL mesh and velocity data couples advection-diffusion-reaction transport to proton NMR detection. The mesh crop is x=`[286.8, 287.5]` and y=`[576.0, 579.0]`, with listed mesh indices inactivated. The source sets exo/endo channel rates to `k1=2.0` and `k2=1.0` (the source comment gives `mol/(L*s)`; the implemented bimolecular terms require `L/(mol*s)` for concentrations in mol/L) and diffusion to `1e-7` without stating its unit. The initial field places 0.50 and 0.25 in reactant components at cells 1240 and 1246; their units are not given. Chemistry uses 501 steps of 20 seconds, while cellwise concentrations are interpolated with `makima` for spin evolution.

The spin system uses `sys.magnet=14.1` (the file does not annotate a field unit) and a rectangular coil-selection mask defined by `287.0<x<287.3` and `577.0<y<577.5`. Proton acquisition parameters are offset 2328, sweep 3500, and 1024 steps; the file does not state units for offset or sweep, and sets the sampling interval to `1/sweep`. It applies a `pi/2` excitation and advances the transport-, relaxation-, and concentration-dependent spin evolution with `step` using left/right interval generators.

## Observable

Acquisitions begin at every 25th chemistry-grid index (21 sampled starts, from 0 through 10,000 seconds). Each FID is apodised, zero-filled to 16384 points, Fourier transformed, and displayed as a real-intensity waterfall against time (seconds) and chemical shift (ppm); intensity is labelled a.u. The source header estimates days and says GPU execution is much faster; this is not a benchmark. The script defaults to `sys.enable={'greedy'}`; GPU transfers are conditional on `gpu` being enabled.

## Source

https://github.com/IlyaKuprov/Spinach/blob/main/examples/microfluidics/reacting_flow_nmr.m
