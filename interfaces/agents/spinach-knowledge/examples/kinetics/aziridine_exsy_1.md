# examples/kinetics/aziridine_exsy_1.m

## Purpose and callable context

The no-argument MATLAB example simulates a two-dimensional NOESY/EXSY experiment on phenylaziridine in a relatively slow-exchange regime. Its source describes scalar relaxation of the second kind (SRSK) from 14N, with first-kind scalar relaxation not manifesting in this case, and says it is set to reproduce Figures 1b and 4b of the cited paper. The source estimates a calculation time of minutes. These are source comments, not a validation of the plotted simulation against the paper.

## Spin and exchange model

The system contains two 10-spin conformational blocks (20 spins total), each with a 14N at spin positions 4 and 14 and otherwise 1H nuclei. The coordinates are explicitly labelled Angstrom in the source; the Zeeman and coupling tensors, including the 14N quadrupolar tensors, are entered as matrices. The source says parameters other than isotropic chemical shifts, exchange rates and correlation times come from a DFT calculation. Isotropic shifts are set in an explicit two-conformer list and the two blocks exchange through `inter.chem.parts`.

The kinetic inputs are `kplus=4` and `kminus=20`, with rate matrix `[-kplus kminus; kplus -kminus]` and concentration weights `[kminus kplus]`. The source does not annotate units for these rate values. Relaxation is configured as Redfield plus SRSK, with SRSK sources at spins 4 and 14, `tau_c={50e-12 50e-12}`, zero equilibrium, and `rlx_keep='kite'`. The spin basis uses `sphten-liouv`, `IK-1`, scalar-coupling connectivity, inter-level 4 and proximity level 3; Krylov is disabled, greedy is enabled, and the proximity cutoff is 10.0.

## Sequence, observable, and limits

The concentration-aware initial state is 1H longitudinal magnetisation. The example calls `liquid` with `@noesy` in NMR mode, uses a mixing-time parameter of 0.800, 128 points and 512 zero-filled points in each dimension, and sets the displayed axes to ppm. The numeric offset and sweep are 2000 and `[4000 4000]`; their units are not specified in the source. It applies squared-cosine apodisation to both cosine and sine channels, forms the States signal, Fourier-transforms both dimensions, and passes the negative real spectrum to `plot_2d`. No numerical spectrum, cross-peak intensities, or independent reproduction result is claimed here.

## Reference

[Article DOI](https://doi.org/10.1002/ange.201410271)

[MATLAB source](https://github.com/IlyaKuprov/Spinach/blob/main/examples/kinetics/aziridine_exsy_1.m)
