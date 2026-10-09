# examples/microfluidics/plain_nmr.m

Source: [examples/microfluidics/plain_nmr.m](https://github.com/IlyaKuprov/Spinach/blob/main/examples/microfluidics/plain_nmr.m)

## System and model

This is a homogeneous liquid-state proton NMR calculation for a Diels–Alder reaction mixture assembled by `dac_reaction()`. The source sets the five concentration entries to `[1 1 1 1 0]` and comments that the four chemical species have equal concentrations with no solvent. Despite its location alongside microfluidics examples, this script includes neither spatial transport nor chemical kinetics. It explicitly removes the builder's reaction records before calling the linear acquisition context.

The field parameter is `sys.magnet=14.1`; this file does not annotate its unit. Greedy parallelisation is enabled. Preparation uses concentration-weighted `state`; detection uses unweighted `coil_state` for `1H`. The solvent block is spin-free and has zero concentration, so it contributes neither preparation nor signal.

## Acquisition and processing

The single-dimensional acquisition uses offset `2328`, sweep `3500`, `4096` points and zero filling to `16384`; the plotted axis unit is ppm and the axis is inverted. The source leaves decoupling empty. `liquid(spin_system,@acquire,parameters,'nmr')` generates the FID, which is apodised with an exponential parameter of `6`, Fourier transformed with `fftshift(fft(fid,parameters.zerofill))`, and plotted as the real spectrum.

## Scope

This is a single liquid-state spectrum calculation: it does not propagate reaction kinetics or a spatially varying signal.

Zero track elimination is explicitly enabled with `zte` in `sys.enable`.
