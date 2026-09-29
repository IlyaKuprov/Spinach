# examples/nmr_paramag/calbindin/tm_1igv_point_fit.m

- MATLAB implementation: [examples/nmr_paramag/calbindin/tm_1igv_point_fit.m](https://github.com/IlyaKuprov/Spinach/blob/main/examples/nmr_paramag/calbindin/tm_1igv_point_fit.m)

- Signature: `tm_1igv_point_fit()`

This example calls a point-electron inverse PCS model to estimate an electron location and susceptibility tensor from the experimental PCS values in `tm_1igv_pcs.mat`. The source comments credit the experimental dataset to Gottfried Otting of the Australian National University.

The function loads `expt_pcs` and coordinate arrays `x`, `y`, and `z`, then calls `ippcs([x y z],[-5 5 -15],expt_pcs)`. The returned values are `mxyz`, `chi`, and `pred_pcs`; the three-number argument is the initial location supplied to the solver. The source does not state coordinate units or a nuclear-isotope assignment, so neither is specified here.

The diagnostic plot places experimental PCS on the x axis and predicted PCS on the y axis, both labelled in ppm. It uses blue circle markers, a red identity line, and shared axis limits spanning the combined experimental and predicted values. The function prints the susceptibility tensor and point-electron location. It does not save a result file or report a fit statistic in this script.
