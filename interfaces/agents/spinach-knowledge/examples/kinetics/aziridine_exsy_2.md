# examples/kinetics/aziridine_exsy_2.m

## Purpose and callable context

The no-argument MATLAB example simulates a two-dimensional NOESY/EXSY experiment on phenylaziridine in the intermediate-exchange regime. The source states that lines broaden, SRSK is induced by 14N, and first-kind scalar relaxation (SRFK) must also be included. It is set to reproduce Figures 1a and 4a of the cited paper; the source estimates hours of calculation time, faster on a GPU. Those comments describe the intended example, not a validation of the simulation output.

## Spin and exchange model

The 20-spin model consists of two 10-spin conformational blocks with 14N at spin positions 4 and 14 and 1H at the other positions. Coordinates are labelled Angstrom, and Zeeman and coupling tensors (including the 14N quadrupolar tensors) are provided as matrices. As the source notes, all parameters except isotropic shifts, exchange rates and correlation times come from a DFT calculation. Shifts are supplied as two ten-spin lists; `inter.chem.parts` couples the two blocks.

Both kinetic rates are `1.2e3`; the source uses `[-kplus kminus; kplus -kminus]` and concentration weights `[kminus kplus]`. No units are attached to these rate values in the source. Relaxation is configured as Redfield, SRFK and SRSK; SRSK sources are spins 4 and 14, `tau_c={50e-12 50e-12}`, and SRFK uses `srfk_tau_c={[1.0 1/kplus]}` with selected `srfk_mdepth` entries. The basis is `sphten-liouv` / `IK-1`, scalar-coupling connectivity, inter-level 4 and proximity level 3; Krylov is disabled, greedy is enabled, and the proximity cutoff is 10.0.

## Sequence, observable, and limits

A concentration-aware 1H longitudinal-magnetisation state is propagated by `liquid` with `@noesy` in NMR mode. The mixing-time parameter is 0.800; acquisition uses 128 points and 512 zero-filled points per dimension, with axes requested in ppm. Offset and sweep are 2000 and `[4000 4000]`; the source does not state their units. The cosine and sine channels receive squared-cosine apodisation, are combined as a States signal, and are Fourier-transformed in both dimensions. The negative real spectrum is plotted by `plot_2d`. The page does not report numerical peak positions/intensities or claim that the target paper figures were independently reproduced.

## Reference

[Article DOI](https://doi.org/10.1002/ange.201410271)

[MATLAB source](https://github.com/IlyaKuprov/Spinach/blob/main/examples/kinetics/aziridine_exsy_2.m)
