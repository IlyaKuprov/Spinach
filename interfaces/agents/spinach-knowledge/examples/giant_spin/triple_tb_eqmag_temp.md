# examples/giant_spin/triple_tb_eqmag_temp.m

- MATLAB implementation: [examples/giant_spin/triple_tb_eqmag_temp.m](https://github.com/IlyaKuprov/Spinach/blob/main/examples/giant_spin/triple_tb_eqmag_temp.m)

## Purpose

The function `triple_tb_eqmag_temp()` calculates temperature-dependent equilibrium magnetisation and the corresponding theoretical chi-T curve for a triangular three-terbium complex. The source links the example to Figures S27 and S28 in the Supplementary Information of [the cited paper](https://doi.org/10.1002/chem.201703842). The J=6 ground-term ligand-field parameters and g-tensor are attributed to SINGLE_ANISO in MOLCAS; the source estimates hours of calculation time.

## Physical model and parameters

The system has three `E13` sites representing J=6 Tb centres. As in the companion field scan, the source builds each g tensor from a site eigenvalue triplet and orientation matrix via `g=U*diag(g_values)*U'`. The triplets are `[1.497075749,1.495252923,1.481370349]`, `[1.496540374,1.494559188,1.482686210]`, and `[1.497265940,1.494866241,1.481858735]`. Orientations, Tb coordinates, and site-specific rank-2, rank-4 and rank-6 Stevens coefficients are defined in the source; it does not state units for the g eigenvalues or coordinates. Exchange is set as `J=icm2hz(0.003)`, with spin-orbit corrections to dipolar couplings enabled by `sys.enable={'sodd'}`. Stevens coefficients are converted with `icm2hz` and `stev2sph`, then rotated by Wigner matrices and passed through `inter.giant.coeff`.

## Calculation and output

The calculation uses the `zeeman-hilb` basis without approximation and the `leb_2ang_rank_11` powder grid. It fixes `sys.magnet=0.1` (the source does not label a unit for this setting) and scans temperature in Kelvin: 1-10 in steps of 1, 15-45 in steps of 5, 50-200 in steps of 10, then 250 and 300. At each temperature it updates `inter.temperature`, calls `eqmag`, and stores the Z component `Mz=mag(3)`. The source computes `chi_theo=0.5585*T.*(Mz/sys.magnet)`, labelled in cm^3 K/mol, and plots that curve with experimental points loaded as `temperature` and `chiT` from `triple_tb_eqmag.mat`.

The plotted experimental comparison depends on that external MAT file. The source defines the calculation and plot but provides no numerical result in the page itself and does not save a result file.
