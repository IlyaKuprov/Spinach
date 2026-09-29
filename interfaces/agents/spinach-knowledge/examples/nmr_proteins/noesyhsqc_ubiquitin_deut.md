# examples/nmr_proteins/noesyhsqc_ubiquitin_deut.m

[MATLAB source](https://github.com/IlyaKuprov/Spinach/blob/main/examples/nmr_proteins/noesyhsqc_ubiquitin_deut.m)

This is a simulated 1H-1H-15N NOESY-HSQC spectrum of 15N-labelled ubiquitin. The source states 900 MHz and a 90 ms mixing time, assumes no 13C labelling, and represents deuterium nuclei explicitly as spin-1 particles. Its stated cost is a week on 32 cores and 512 GB of RAM.

## Protein model and sequence settings

The model uses all atoms from 1D3Z.pdb / 1D3Z.bmrb. It deuterates the listed positions HA, HB, HB1, HB2, HB3, HG, HG1, HG2, HG3, HD, HD1, HD2, HD3, HE, HE1, HE2, HE3, HZ, HZ1, HZ2, HZ3, HH, HH1, HH2, and HH3; the code then removes 13C spins. The field parameter is sys.magnet=21.1356, without a unit annotation in the code. Interaction/proximity cutoffs are 2.0/4.0. Relaxation settings are Redfield plus T1/T2, r1_rates set to 100 for 2H and 0 otherwise, rlx_keep='kite', zero equilibrium, and tau_c=1e-8; the example does not annotate units for these numeric relaxation parameters. The basis is IK-1 sphten-liouv with scalar-coupling connectivity and inter/proximity levels 4/3; the code disables asyredf.

The call uses liquid with noesyhsqc, spin order 1H, 15N, 1H, tmix=0.090 (90 ms in the source comment), and J=90.0. It sets point counts [128 64 128], zero-fill sizes [512 256 512], offsets [4250 -10600 4250], and sweep settings [10750 3000 10750]; units are not annotated for J, offsets, or sweeps, while the axes are specified in ppm. The example gives no RF pulse widths or phases, contact time, or rotor parameters.

## Simulation and output

Four phase-cycle FIDs are squared-cosine apodised. The code combines conjugate components to form absorption parts through the F3 and F2 transforms, then Fourier-transforms F1 and plots the negative real 3D spectrum. The PDB/BMRB files initialise the simulated model; no experimental spectrum is input and no measured coherence data or DOI citation is reported. The example identifies the 1H-1H-15N sequence dimensions and spin-1 deuterons, but does not label explicit coherence orders.
