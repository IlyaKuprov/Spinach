# examples/nmr_paramag/calbindin/tm_1igv_distr_fit.m

- MATLAB implementation: [examples/nmr_paramag/calbindin/tm_1igv_distr_fit.m](https://github.com/IlyaKuprov/Spinach/blob/main/examples/nmr_paramag/calbindin/tm_1igv_distr_fit.m)

- Signature: `tm_1igv_distr_fit()`

## Purpose and data

This inverse-problem example estimates a spatial distribution of unpaired-electron density from experimental pseudocontact shifts (PCS). The source credits the experimental data to Gottfried Otting (Australian National University). It reads the processed 1IGV PDB, the PCS positions and measurements from `tm_1igv_pcs.mat` (variables `x`, `y`, `z`, and `expt_pcs`), and an effective susceptibility tensor `chi` from `tm_1igv_chi_eff.mat`. The PDB atom coordinates are passed as structural geometry; the source does not state a coordinate unit or a nuclear-isotope list.

## Density reconstruction workflow

The `ipcs` model is configured with equation `kuprov`, a plot selection of diagnostics, density, molecule, tightzoom, and box, box centre `[3.5 17.0 16.1]`, and box size `[7.0 7.0 7.0]`. It passes all PDB atom coordinates as `xyz_all`, sets six margins to 50, confines the reconstruction with `[1.0 3.0]`, uses sharpening value `2e3`, and supplies the measured PCS, their coordinates, and the loaded `chi`. GPU execution is enabled.

The script refines at grid sizes 64, 128, 256, and 384, calling `ipcs(parameters,n,0.34)` at each size and passing each returned density cube forward as the next initial guess. After the final pass, it calls `chi_eff(source_cube,ranges,[x y z],expt_pcs)`, displays the resulting effective susceptibility tensor, and saves it back to `tm_1igv_chi_eff.mat` as `chi`.

## Scope

The file specifies a reconstruction procedure and saved tensor update; it does not report a numerical fit score or establish agreement with the experimental PCS. The units and sign convention for the loaded tensor, coordinates, and PCS are not documented in this script and are therefore left unstated here.
