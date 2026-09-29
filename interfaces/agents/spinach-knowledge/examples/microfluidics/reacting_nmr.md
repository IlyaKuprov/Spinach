# examples/microfluidics/reacting_nmr.m

## Model and kinetics

This is a homogeneous, non-spatial reaction-NMR calculation: it integrates two competing cycloaddition channels, then couples the interpolated concentrations to chemical reaction generators during repeated proton pulse-acquire acquisitions. There is no flow or diffusion. The source comments give `k1=0.5` toward exo and `k2=0.1` toward endo, with bimolecular rate units `L/(mol*s)` for mol/L concentrations and seconds; the initial concentrations are `[0.6; 0.5; 0; 0; 18.1] mol/L`. Concentration kinetics run for 20 seconds in 200 `LG4` steps; makima interpolation supplies concentrations between those points. The file header identifies Redfield relaxation and describes runtime as hours, much faster on GPU; those are source descriptions, not measurements here.

## NMR acquisition and output

At `sys.magnet=14.1`, the script models `1H` with offset 2328, sweep 3500, and 4096 evolution steps per acquisition; units for these parameter values are not stated in the file. It uses relaxation, a `pi/2` pulse, and reaction-dependent left/right interval generators in a two-point Lie quadrature `step`. Acquisitions start at each integer second from 0 to 18. The 19 FIDs are apodised, zero-filled to 16384 points, and Fourier transformed. The displayed waterfall uses real intensity (a.u.), chemical shift (ppm), and time (seconds); the separate concentration plot shows the first four components in mol/L and omits solvent.

## Source

https://github.com/IlyaKuprov/Spinach/blob/main/examples/microfluidics/reacting_nmr.m
