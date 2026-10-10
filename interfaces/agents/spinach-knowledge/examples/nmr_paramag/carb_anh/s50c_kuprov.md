# examples/nmr_paramag/carb_anh/s50c_kuprov.m

- Signature: `s50c_kuprov()`
- Source: [s50c_kuprov.m](https://github.com/IlyaKuprov/Spinach/blob/main/examples/nmr_paramag/carb_anh/s50c_kuprov.m)

## S50C distributed-density reconstruction

This example treats the human carbonic anhydrase II S50C PCS data as a distributed electron-density inverse problem, rather than a point-centre fit. The source cites the [study](https://doi.org/10.1039/c6sc03736d) and the [PCS analysis tutorial](https://spindynamics.org/wiki/index.php?title=Pseudocontact_shift_analysis).

It loads `expt_pcs`, `xyz`, and `xyz_all` from `s50c_expt.mat`, plus an effective susceptibility tensor from `s50c_chi_eff.mat`. The `kuprov` equation is solved by `ipcs` on grids of 64, 128, and 256 points per dimension, using the preceding density cube as the next guess and passing `0.23` as the third solver argument. The configured reconstruction region is centred at `[-27.4 13.3 18.8]`, has size `[25.0 25.0 25.0]`, and uses margins of 50, confinement `[3.0 12.0]`, and sharpening 1.0. Density/molecule/diagnostic plots are requested and GPU execution is enabled. Finally, `chi_eff` evaluates and displays the tensor for the final cube.

The script does not name the observed nucleus or state field, temperature, coordinate units, or tensor units; the numerical region parameters therefore have no source-stated physical units here. It produces a reconstruction and tensor display, not a spectrum, and does not explicitly save the final cube or tensor.
