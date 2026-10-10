# examples/nmr_nucleic/rna_noesy_theo.m

- Signature: `rna_noesy_theo()`
- Source: [examples/nmr_nucleic/rna_noesy_theo.m](https://github.com/IlyaKuprov/Spinach/blob/main/examples/nmr_nucleic/rna_noesy_theo.m)

## What this example models

This wrapper simulates a theoretical 1H-1H NOESY spectrum for the example RNA provided by Gerhard Wagner's group at Harvard University. It imports `example.pdb` and `example.txt` with `nuclacid`, deletes shifts, and lists exchangeable protons for deuteration: GUA:H1, GUA:H21, GUA:H22, CYT:H41, CYT:H42, URI:H3, ADE:H61, and ADE:H62. After system construction it removes 13C and 15N spins, reflecting the source comment's assumption that the RNA is unlabelled; the simulated spin set is therefore proton-only.

The source sets the field to 17.62 T. Relaxation is Redfield, with `rlx_keep='kite'`, zero equilibrium, and `tau_c={3e-9}`. The basis uses `sphten-liouv`, `IK-1`, scalar-coupling connectivity, interaction level 5, and proximity level 3. Interaction/proximity cutoffs are 1.0 and 5.0. Krylov propagation and colorbar are disabled; propagator caching and greedy are enabled. Units for these numeric settings, including `tau_c`, are not stated in this wrapper. Zero track elimination is also explicitly enabled with `zte` in `sys.enable`.

## Sequence parameters and processing

The mixing parameter is `tmix=0.200`; the wrapper does not annotate its unit. It sets offset 3473, sweeps `[7500 7500]`, acquired points `[512 1024]`, zero-fill sizes `[1024 4096]`, a 1H spin label, and ppm axis labels. As in the source, the offset and sweep values are given without explicit units. The initial density operator is built from proton `Lz` using the `cheap` state option.

The simulated signal is returned by `liquid(spin_system,@noesy,parameters,'nmr')`. The cosine and sine components are squared-cosine apodised, transformed along F2, combined as `f1_cos-1i*f1_sin` (the States signal in the source comment), then transformed along F1. The wrapper plots the negative real part of that spectrum.

## Scope and limits

The wrapper calls the sequence through `@noesy` but does not provide its internal pulse, gradient, or receiver-phase details; those are not inferred. Its comment estimates hours of calculation; that is not a measured runtime. This is a simulation driven by supplied RNA structure/assignment files, not a readout of an experimental spectrum. Source comments credit Shunsuke Imai, Scott Robson, Gerhard Wagner, Zenawi Welderufael, and Ilya Kuprov.
