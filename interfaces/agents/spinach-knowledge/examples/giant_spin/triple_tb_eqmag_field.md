# examples/giant_spin/triple_tb_eqmag_field.m

- MATLAB implementation: [examples/giant_spin/triple_tb_eqmag_field.m](https://github.com/IlyaKuprov/Spinach/blob/main/examples/giant_spin/triple_tb_eqmag_field.m)

## Purpose

The function `triple_tb_eqmag_field()` calculates the field dependence of the equilibrium magnetisation of a triangular three-terbium complex. The source associates the calculation with Figures S27 and S28 in the Supplementary Information of [the cited paper](https://doi.org/10.1002/chem.201703842). It says the J=6 ground-term ligand-field parameters and g-tensor were computed with SINGLE_ANISO in MOLCAS; estimated calculation time is hours.

## Physical model and parameters

The three giant-spin sites are entered as `E13` isotopes, each representing a J=6 Tb centre. Site g tensors are assembled as `g=U*diag(g_values)*U'`. The three eigenvalue triplets in the source are `[1.497075749,1.495252923,1.481370349]`, `[1.496540374,1.494559188,1.482686210]`, and `[1.497265940,1.494866241,1.481858735]`; the corresponding orientation matrices and coordinates are defined in the source. The source does not assign units to those tensor or coordinate entries.

Equal pairwise scalar exchange couplings are set with `J=icm2hz(0.003)`; the stated conversion uses the Spinach NMR convention. `sys.enable={'sodd'}` enables the source-commented spin-orbit corrections to dipolar couplings. Site-specific Stevens coefficients of ranks 2, 4 and 6 are converted with `icm2hz`, converted to irreducible spherical tensors with `stev2sph`, and rotated with Wigner matrices before assignment to `inter.giant.coeff`. The source also defines site Euler data in `inter.giant.euler`.

## Calculation and output

The basis is `zeeman-hilb` with `bas.approximation='none'`. At fixed `inter.temperature=2.0` K, the script evaluates `eqmag` over the field array `[0.01,0.1,0.2,...,1.5,2,3,4,5,6]`; the plot labels field in Tesla. For each field it creates the spin system and basis, takes `mag(3)` as the Z component, and plots magnetisation labelled in Bohr magnetons. The powder grid is `leb_2ang_rank_11`. The theory curve is overlaid with points loaded from `triple_tb_eqmag.mat` (`field` and `magn`).

This is a script-defined calculation and plotting workflow, not a reported numerical result. The comparison requires the external MAT file; the source does not save a result file.
