# examples/nmr_proteins/hsqc_ubiquitin_b.m

[MATLAB source](https://github.com/IlyaKuprov/Spinach/blob/main/examples/nmr_proteins/hsqc_ubiquitin_b.m)

This example simulates a 1H-15N HSQC of human ubiquitin. Its source comment specifies no 1H decoupling in F1 and no 15N decoupling in F2, retaining nitrogen-proton multiplicity in both dimensions; the code configures 13C decoupling in both dimensions. The source credits Zenawi Welderufael, Luke Edwards, and Ilya Kuprov.

## Protein model and sequence settings

The model is imported from 1D3Z.pdb and 1D3Z.bmrb with the backbone-HSQC selection and the noshift='delete' option. The field parameter is sys.magnet=14.1; the example does not annotate its unit. Interaction and proximity cutoffs are 5.0 and 4.0. It uses an IK-1 sphten-liouv basis with scalar-coupling connectivity, inter level 4 and proximity level 1.

The simulation calls liquid with the hsqc sequence and spins 15N and 1H, J=90, sweep settings [2400 4800], offsets [-7000 4500], 256 points per dimension, and zero-fill sizes [1024 1024]. The code does not state units for J, sweep, offsets, or magnet setting; it does state ppm for plotted axes. No RF pulse widths or phases, contact time, or rotor parameters are specified.

## Simulation and output

The generated simulated FIDs are squared-cosine apodised. The code Fourier-transforms F2, combines the two components as f1_pos+conj(f1_neg) for the States signal, then transforms F1 and plots the magnitude spectrum. The PDB/BMRB files initialise the simulated protein model; the example does not load an experimental spectrum or report measured coherence data. No DOI citation appears in the source.
