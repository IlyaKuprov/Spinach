# examples/fitting/nmr_rdc/saupe_example.m

- MATLAB implementation: [examples/fitting/nmr_rdc/saupe_example.m](https://github.com/IlyaKuprov/Spinach/blob/main/examples/fitting/nmr_rdc/saupe_example.m)

- Signature: `saupe_example()`

## Purpose

This example fits a Saupe order matrix from residual dipolar coupling (RDC) measurements and compares the measured couplings with values back-calculated from that matrix. The source header describes the measurements as NH RDC data and credits Andras Boeszoermenyi, Thibault Viennet, and Hari Arthanari.

## Inputs and fit parameterisation

Run `saupe_example()` from a working directory containing `protein.pdb` and `rdc_data.mat`. The function has no input arguments or declared return values. It reads the structure with `read_pdb_pro('protein.pdb',1)`; from the MAT file it loads `aa_num`, `rdc_hz`, `isotope_a`, `isotope_b`, `atom_a`, and `atom_b`. For each measured coupling, it selects the two coordinates matching the residue number and atom labels, and passes the isotope pairs, coordinate pairs, and measured RDC vector to `rdc_fit(isotopes,xyz,rdc_hz)`. Thus the example's explicit fit inputs are the observed couplings and their isotope/geometry assignments; it does not specify fit initialisation, objective details, or restraints in this script.

## Back-calculation and output

The fitted matrix `S` is printed. Each coupling is then recomputed by `xyz2rdc` using its isotope pair and coordinates with the model selector `{S,'saupe'}`. The figure plots experimental RDC (horizontal axis) against theoretical RDC (vertical axis), both labelled in Hz, with red data points, a blue identity line, square axes, and fixed limits from −30 to 30 Hz. The script does not report a fit score or save the figure or matrix.

The atom-selection masks are used directly to index the PDB outputs; the example does not check that each mask selects exactly one atom. The result therefore depends on the supplied structure labels and RDC table agreeing. The plotted range is fixed by the script and should not be read as a stated range for the input measurements.
