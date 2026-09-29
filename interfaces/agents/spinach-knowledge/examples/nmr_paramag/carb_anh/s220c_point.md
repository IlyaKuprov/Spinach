# examples/nmr_paramag/carb_anh/s220c_point.m

- Signature: `s220c_point()`
- Source: [s220c_point.m](https://github.com/IlyaKuprov/Spinach/blob/main/examples/nmr_paramag/carb_anh/s220c_point.m)

## S220C point-centre fit

This human carbonic anhydrase II example fits the S220C PCS data with a point paramagnetic-centre model. Its source cites the [study](https://doi.org/10.1039/c6sc03736d) and the [PCS analysis tutorial](https://spindynamics.org/wiki/index.php?title=Pseudocontact_shift_analysis).

The script loads experimental PCS values and coordinates from `s220c_expt.mat`, then calls `ippcs` with the coordinate data and the supplied vector `[-14 -26 4]`. The fit returns the point location `mxyz`, susceptibility tensor `chi`, and predicted PCS values. It plots predicted against experimental PCS in ppm with a diagonal reference and displays `chi` and `mxyz`; it does not save those results.

The source does not identify the observed nucleus, give field or temperature values, or state units for the coordinate vector, fitted centre, or susceptibility tensor. It describes a PCS fit, not a spectral simulation; those unspecified quantities should not be inferred from this script.
