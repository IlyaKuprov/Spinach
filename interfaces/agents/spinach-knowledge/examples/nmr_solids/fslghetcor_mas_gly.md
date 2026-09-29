# examples/nmr_solids/fslghetcor_mas_gly.m

- Signature: fslghetcor_mas_gly()

## Purpose and spin model

This script sets up an FSLG-HETCOR calculation for alpha-glycine powder under MAS. It imports a PCM-DFT spin system from the relative-path glycine.log file with two ¹³C and five ¹H labels, then assigns isotropic chemical-shift values: 176.4 and 43.6 for the two carbons, 2.6 and 3.8 for the alpha protons, and 8.0 for each of the three nitrogen-bound protons. The source labels these as alpha-glycine shifts but does not annotate their units. The remaining spin-system interactions come from the parsed input log; this script does not print their tensor values. It sets sys.magnet to 9.4 and comments that this is a 400 MHz spectrometer. Interactions below 200 Hz are excluded using a 2π×200 cutoff. The basis is sphten-liouv with IK-0 approximation and inter-level 4.

Initial state is ¹H Lz and detection is on ¹³C L+. The experiment invokes fslghetcor through singlerot, so the method is heteronuclear correlation with frequency-switched Lee–Goldburg (FSLG) rather than an HMQC sequence.

## MAS, transfer, and acquisition settings

The code sets the MAS-rate value to 10000 (without a unit comment), uses the rep_2ang_100pts_sph powder grid and maximum rank 7, and passes high-power RF amplitude 83e3 Hz, CP powers [60e3, 50e3] Hz, and a contact duration of 1e-4 seconds. It uses four FSLG blocks. The two CP amplitudes are code-set parameters; the script does not report a separately measured or calibrated Hartmann–Hahn match.

The acquisition array is ordered as F1 then F2: 128 and 512 points, zero-filled to 512 and 2048. The F2 sweep is 1/(33e-6); the F1 sweep is calculated from the high-power amplitude and block count with the explicit magic-angle scaling in the source. The F1 offset is set to zero after initialisation, while the F2 offset remains 1e4 Hz.

## Processing and plotted observable

The simulation returns cosine and sine channels, each apodised with squared-cosine windows in both dimensions. The script Fourier-transforms the indirect dimension for both channels, forms States quadrature from their real parts, then Fourier-transforms the direct dimension. It plots the real 2D spectrum. The source comment estimates hours on an NVIDIA Tesla A100 and longer on CPU, but GPU arithmetic is commented out in this function. The wrapper call defines the FSLG-HETCOR simulation settings but does not document unlisted sequence internals. The source reports no measured spectrum or transfer-efficiency result. No DOI is given in the source.

Source: https://github.com/IlyaKuprov/Spinach/blob/main/examples/nmr_solids/fslghetcor_mas_gly.m
