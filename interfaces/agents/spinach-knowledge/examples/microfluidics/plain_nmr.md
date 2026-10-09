# examples/microfluidics/plain_nmr.m

Source: [examples/microfluidics/plain_nmr.m](https://github.com/IlyaKuprov/Spinach/blob/main/examples/microfluidics/plain_nmr.m)

## System and model

This is a homogeneous liquid-state proton NMR calculation for a Diels–Alder reaction mixture assembled by `dac_reaction()`. The source sets the five concentration entries to `[1 1 1 1 0]` and comments that the four chemical species have equal concentrations with no solvent. Despite its location alongside microfluidics examples, this script includes neither spatial transport nor chemical kinetics. It explicitly removes the builder's reaction records before calling the linear acquisition context.

The field parameter is `sys.magnet=14.1`; this file does not annotate its unit. Greedy parallelisation is enabled. Preparation uses concentration-weighted `state`; detection uses unweighted `coil_state` for `1H`. The solvent retains three protons with T1/T2 relaxation but has zero concentration, so it contributes neither preparation nor signal.

## Acquisition and processing

The single-dimensional acquisition uses offset `2328`, sweep `3500`, `4096` points and zero filling to `16384`; the plotted axis unit is ppm and the axis is inverted. The source leaves decoupling empty. `liquid(spin_system,@acquire,parameters,'nmr')` generates the FID, which is apodised with an exponential parameter of `6`, Fourier transformed with `fftshift(fft(fid,parameters.zerofill))`, and plotted as the real spectrum.

## Scope

This is a single liquid-state spectrum calculation: it does not propagate reaction kinetics or a spatially varying signal.

Zero track elimination is explicitly enabled with `zte` in `sys.enable`.

## Numerical batching

Direct-sum reduction propagates the small acrylonitrile block separately, using the small-matrix `expm` path. The former global reduction batched it with cyclopentadiene and used the rounded Taylor/scaling-and-squaring path. At the default `prop_chop=1e-10`, this produces a relative full-FID difference of about `1.63e-6` despite identical represented descriptors, relaxation, preparation, and detection, and Hamiltonians equal to round-off. Reconstructing the acrylonitrile propagator with the former batch scaling accounts for that difference to `7.4e-14` relative full-FID residual. An independent dense MATLAB `expm` reference on all signal-bearing invariant generator blocks gives whole-FID relative errors of `4.8451e-6` for the former stock trajectory and `4.5633e-6` for the direct-sum trajectory. Thus the migrated trajectory is closer to the reference, not a loss of accuracy. Both trajectories share the larger-block rounded-propagator error; separating acrylonitrile removes its additional stock batched-Taylor error. The two-half-step `expm` construction gives the same reference to machine precision. This is a documented numerical PASS-WITH-REASON, not a literal `1e-10` WP0 match or a changed chemical model; physical parameters and numerical defaults are retained.

Detection and reference operator vectors explicitly use the `exact` method of the four-argument `coil_state` primitive.
