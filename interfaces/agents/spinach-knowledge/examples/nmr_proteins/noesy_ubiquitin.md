# examples/nmr_proteins/noesy_ubiquitin.m

[MATLAB source](https://github.com/IlyaKuprov/Spinach/blob/main/examples/nmr_proteins/noesy_ubiquitin.m)

This example simulates a 1H-1H NOESY spectrum of ubiquitin with a 65 ms mixing time. The source assumes the protein is neither 13C- nor 15N-labelled and estimates hours of calculation time.

## Protein model and sequence settings

The example imports all atoms from 1D3Z.pdb / 1D3Z.bmrb, then removes 13C and 15N spins. It sets sys.magnet=21.1356 (unit not annotated in this source), interaction/proximity cutoffs of 2.0/4.0, and Redfield relaxation with rlx_keep='kite', zero equilibrium, and tau_c=5e-9 (unit not specified). The basis is IK-1 sphten-liouv with scalar-coupling connectivity and inter/proximity levels 4/3.

The sequence uses liquid with noesy; the source sets tmix=0.065, matching its stated 65 ms mixing time. It also sets offset 4250, sweep settings [11750 11750], 512 points per dimension, and zero-fill sizes [2048 2048]. Those numeric settings are not assigned units in the file; the plotted axes are in ppm. RF pulse widths or phases, contact time, and rotor parameters are not specified.

## Simulation and output

The simulated cosine and sine FIDs are squared-cosine apodised and Fourier-transformed in F2. The States signal is formed as f1_cos-1i*f1_sin before the F1 transform, and the plot shows the negative real part of the 2D spectrum. The PDB/BMRB inputs define the simulated protein model; the example does not read an experimental spectrum or report measured coherence data. No DOI citation appears in the source.
