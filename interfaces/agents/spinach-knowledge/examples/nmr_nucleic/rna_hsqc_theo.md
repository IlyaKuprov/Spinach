# examples/nmr_nucleic/rna_hsqc_theo.m

- Signature: `rna_hsqc_theo()`

## Purpose

Simulates a 1H–13C HSQC spectrum for the example RNA molecule supplied by the Wagner group. The source estimates a calculation time of minutes and credits Shunsuke Imai, Scott Robson, Gerhard Wagner, Zenawi Welderufael, and Ilya Kuprov.

## Physical / mathematical content

The RNA spin system is imported from `example.pdb` and `example.txt` with shift deletion enabled. At 17.62 T, the calculation uses an sphten-liouv basis with IK-1 approximation and scalar-coupling connectivity. The selected spins are 13C and 1H; the sequence parameters set J=90, dimension sweeps [2500, 2000], offsets [22000, 6000], and acquisition sizes [128, 256]. Relaxation is damped with diagonal retention, zero equilibrium, and damping rate 5.0.

## Numerical / algorithmic content

The source disables Krylov propagation, enables the greedy option, and builds the basis through `create` and `basis`. It simulates `@hsqc`, applies cosine apodisation to the positive- and negative-frequency signals, Fourier-transforms F2, combines the components as a States signal, and Fourier-transforms F1 before plotting the real spectrum. Both dimensions are zero-filled to 1024 points.

## Implementation structure

The workflow imports the RNA, configures system tolerances and basis, sets sequence and relaxation parameters, calls `liquid(spin_system,@hsqc,parameters,'nmr')`, processes the two signal components, and plots the resulting 2D spectrum.
