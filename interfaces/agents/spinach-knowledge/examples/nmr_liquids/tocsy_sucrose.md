# examples/nmr_liquids/tocsy_sucrose.m

- Signature: `tocsy_sucrose()`
- Source: [examples/nmr_liquids/tocsy_sucrose.m](https://github.com/IlyaKuprov/Spinach/blob/main/examples/nmr_liquids/tocsy_sucrose.m)

## What it models

This is a simulated liquid-state, two-dimensional `1H` TOCSY spectrum of sucrose. The source comment says that its magnetic parameters were computed with DFT and gives a calculation-time estimate of seconds; these are source comments, not a reported run or measured spectrum. It parses `../standard_systems/sucrose.log` with `gparse`, then calls `g2spinach` for hydrogen nuclei labelled `1H`. The wrapper does not load an experimental FID or spectrum. It does not establish whether any separately supplied input values were measured; its spin-system input is labelled vacuum DFT in the source.

## Spin system and approximation

- The conversion call supplies the numeric argument `31.8`; the source does not label its units or explain its meaning. The source sets `options.min_j=1.0` and `sys.magnet=5.9`; units are not stated in this wrapper.
- The basis is `sphten-liouv`, with `IK-2` approximation, `scalar_couplings` connectivity, and proximity level `1`.
- It enables `greedy`, disables `krylov`, and sets `sys.tols.prox_cutoff=4.0`. The initial density operator is `state(spin_system,'Lz','1H')`.
- The wrapper assigns no relaxation parameters. Do not infer an additional relaxation model from this file.

## Sequence call and processing

The source sets `parameters.tmix=0.100`, `parameters.lamp=1e4`, `parameters.offset=800`, `parameters.sweep=[1700 1700]`, `parameters.npoints=[512 512]`, and `parameters.zerofill=[2048 2048]`. It sets observed spins to `{'1H'}` and display-axis units to `ppm`. Aside from that explicit axis label, this wrapper does not state units for these numeric parameter values; `lamp` is not defined further here.

The simulation call is `liquid(spin_system,@tocsy,parameters,'nmr')`. This delegates the sequence to the `tocsy` helper. This wrapper does not specify pulse angles/phases, gradients, receiver phases, or internal phase cycling, so those details cannot be described from this source alone.

The returned cosine and sine FIDs are each apodised in both dimensions with `sqcos`. The code Fourier-transforms along the dimension it labels F2, takes the imaginary component from the cosine channel and the real component from the sine channel, combines them as `f1_cos-1i*f1_sin` (States signal), then Fourier-transforms along F1. It plots `abs(spectrum)` with `plot_2d`; the function has no output argument and writes no spectrum file in this wrapper.
