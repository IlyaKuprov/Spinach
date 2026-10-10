# examples/nmr_nucleic/rna_hsqc_theo.m

- Signature: `rna_hsqc_theo()`
- Source: [examples/nmr_nucleic/rna_hsqc_theo.m](https://github.com/IlyaKuprov/Spinach/blob/main/examples/nmr_nucleic/rna_hsqc_theo.m)

## What this example models

This example simulates a theoretical 1H-13C HSQC spectrum for the example RNA supplied by the Wagner group. It imports structural and assignment information from `example.pdb` and `example.txt` through `nuclacid`, with `options.noshift='delete'` and an empty deuteration list. These files define the modeled RNA system; the wrapper does not load an experimental spectrum.

The field is 17.62 T. The basis is `sphten-liouv`, approximation `IK-1`, scalar-coupling connectivity, interaction level 4 and proximity level 1. Interaction and proximity cutoffs are both set to 5.0 (the wrapper does not state units). Krylov propagation and colorbar are disabled, and the greedy option is enabled. Relaxation is damped with diagonal retention and zero equilibrium; `damp_rate=5.0`, with no unit annotation in this file. Zero track elimination is also explicitly enabled with `zte` in `sys.enable`.

## Sequence parameters and processing

The wrapper supplies `J=90`, sweeps `[2500 2000]`, offsets `[22000 6000]`, acquired points `[128 256]`, and zero-fill sizes `[1024 1024]`. It labels the spins `{'13C','1H'}`, sets F1 proton decoupling and F2 carbon decoupling, and labels the axes ppm. The numeric J, sweep, offset, cutoff, and damping-rate values are reproduced as set by the source; this wrapper does not attach units to them.

The simulated signal comes from `liquid(spin_system,@hsqc,parameters,'nmr')`. The returned positive and negative components are cosine-apodised, Fourier transformed along the F2 dimension, combined as `f1_pos+conj(f1_neg)` (the States combination in the source comment), then Fourier transformed along F1. The plotted quantity is the real part of the resulting spectrum.

## Scope and limits

The wrapper specifies the sequence by the `@hsqc` handle but does not expose the pulse program's internal pulse/gradient/phase-cycle details; none are inferred here. The source comment estimates minutes of calculation, which is an estimate rather than a measured runtime. Source comments credit Shunsuke Imai, Scott Robson, Gerhard Wagner, Zenawi Welderufael, and Ilya Kuprov.
