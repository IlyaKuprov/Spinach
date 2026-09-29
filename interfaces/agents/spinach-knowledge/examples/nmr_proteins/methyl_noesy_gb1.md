# examples/nmr_proteins/methyl_noesy_gb1.m

[MATLAB source](https://github.com/IlyaKuprov/Spinach/blob/main/examples/nmr_proteins/methyl_noesy_gb1.m)

This example simulates a 1H-1H NOESY spectrum of GB1 with non-methyl positions deuterated. The source says deuteria remain in the spin system because they belong to the coupling network, while methyl-group rotation is not modelled; it estimates hours of calculation time.

## Protein model and sequence settings

The input is 2N9K.pdb / 2N9K.bmrb with all atoms selected and the non-methyl positions deuterated. The code then removes 13C and 15N spins, stating that the protein is assumed to be unlabelled. The magnet parameter is sys.magnet=21.1356 (unit not annotated in this source). The interaction cutoff is 100, commented as retaining only significant dipolar couplings; the proximity cutoff is 5.0, with the comment “increase till convergence.” The model uses Redfield relaxation, rlx_keep='kite', zero equilibrium, and tau_c=5e-9 (unit not specified), with an IK-1 sphten-liouv scalar-coupling basis at inter/proximity levels 2/2.

The sequence call uses liquid with noesy, tmix=200e-3, offset 750, sweep settings [3000 3000], 512 points per dimension, and zero-fill sizes [2048 2048]. The source does not annotate units for these numeric sequence settings or for tmix; the axes are specified in ppm. It gives no RF pulse widths or phases, contact time, or rotor parameters.

## Simulation and output

The simulation produces cosine and sine FIDs, which are squared-cosine apodised and transformed in F2. Their States combination is f1_cos-1i*f1_sin, followed by the F1 transform; the plotted quantity is the negative real part of the 2D spectrum. The PDB/BMRB files provide the protein-model input, not an experimental spectrum; no experimental coherence observation or DOI citation is supplied in the example.
