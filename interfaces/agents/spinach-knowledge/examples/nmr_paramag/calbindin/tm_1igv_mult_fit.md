# examples/nmr_paramag/calbindin/tm_1igv_mult_fit.m

- Signature: `tm_1igv_mult_fit()`

## Purpose

Recovers a point-electron location and susceptibility tensor from experimental PCS data. The source credits the experimental data to Gottfried Otting (Australian National University).

## Numerical and implementation details

The script loads PCS values and x, y, z coordinates from `tm_1igv_pcs.mat` and calls `ilpcs([x y z],expt_pcs,[0 1 2],[-5 5 -15])` to obtain the fitted location, susceptibility tensor and predicted PCS values. It plots measured versus predicted PCS with a y=x reference line, then displays the recovered tensor and point-electron location.
