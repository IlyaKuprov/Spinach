# examples/nmr_paramag/carb_anh/s220c_kuprov.m

- Function: `s220c_kuprov()`
- Source: [`examples/nmr_paramag/carb_anh/s220c_kuprov.m`](https://github.com/IlyaKuprov/Spinach/blob/main/examples/nmr_paramag/carb_anh/s220c_kuprov.m)

## Purpose

Reconstructs a distributed electron-density model from PCS data for the S220C mutant of human carbonic anhydrase II. The example cites method paper DOI [10.1039/c6sc03736d](https://doi.org/10.1039/c6sc03736d) and the [PCS analysis tutorial](http://spindynamics.org/wiki/index.php?title=Pseudocontact_shift_analysis).

## Inputs and density fit

Loads `expt_pcs`, `xyz`, and `xyz_all` from `s220c_expt.mat`, and an initial effective susceptibility tensor `chi` from `s220c_chi_eff.mat`. It configures `ipcs` with equation `kuprov`, GPU execution, box centre [-16.0, -25.5, 6.0], box size [30.0, 20.0, 25.0], margins 50 in each of six directions, confinement [3.0, 12.0], and sharpening 1.0. The plot requests are diagnostics, density, molecule, tight zoom, and box.

The solver is called successively at grid sizes 64, 128, and 256 with parameter 0.17. Each returned source cube becomes the next call's initial guess. After the final grid, `chi_eff(source_cube,ranges,xyz,expt_pcs)` calculates an effective susceptibility tensor, which the function displays.

## Scope and limitations

The source identifies PCS but does not specify the measured nuclei, field, temperature, coordinate units, or tensor units. The numerical box and confinement settings above are source literals, not unit assignments. The function does not state the fitted density as a numerical result in the file; it displays the updated tensor.
