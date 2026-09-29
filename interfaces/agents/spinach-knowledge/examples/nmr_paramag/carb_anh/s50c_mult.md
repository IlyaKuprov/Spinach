# examples/nmr_paramag/carb_anh/s50c_mult.m

- Signature: `s50c_mult()`
- Source: [s50c_mult.m](https://github.com/IlyaKuprov/Spinach/blob/main/examples/nmr_paramag/carb_anh/s50c_mult.m)

## S50C multipolar PCS fit

This human carbonic anhydrase II S50C example uses a multipolar PCS model with orders 0, 1, and 2; it is distinct from both the S220C point-centre fit and the S50C distributed-density reconstruction. The source cites the [study](https://doi.org/10.1039/c6sc03736d) and the [PCS analysis tutorial](https://spindynamics.org/wiki/index.php?title=Pseudocontact_shift_analysis).

It loads experimental PCS and coordinates from `s50c_expt.mat` and calls `ilpcs` with the order set `[0 1 2]` and supplied centre vector `[-27.0 13.0 18.0]`. The call returns the multipole-centre location `mxyz`, susceptibility tensor `chi`, and predicted PCS values. The example plots predicted versus experimental PCS in ppm against a diagonal reference, then displays the tensor and fitted centre; it does not save these results.

The source does not specify the observed nucleus, field, temperature, or units for the supplied/fitted centre or susceptibility tensor. No spectral simulation is performed by this script.
