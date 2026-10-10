# examples/nmr_paramag/carb_anh/s166c_mult.m

- MATLAB implementation: [examples/nmr_paramag/carb_anh/s166c_mult.m](https://github.com/IlyaKuprov/Spinach/blob/main/examples/nmr_paramag/carb_anh/s166c_mult.m)

- Signature: `s166c_mult()`

The source calls this a multipolar fit for the S166C mutant dataset of human carbonic anhydrase II. It cites https://doi.org/10.1039/c6sc03736d for the system and method, and http://spindynamics.org/wiki/index.php?title=Pseudocontact_shift_analysis for a tutorial. The article gives the literature context; the workflow below describes the script's simulation.

The function loads `expt_pcs` and `xyz` from `s166c_expt.mat`, then calls `ilpcs(xyz,expt_pcs,[0 1 2],[-15.0 -3.0 -10.0])`. The source labels this as the inverse problem and receives `mxyz`, `chi`, and `pred_pcs` from the call (with the third output discarded); the supplied order selector is [0, 1, 2], and the final three-number argument is passed as the initial location. The source does not identify coordinate units or nuclear isotopes.

The plot compares experimental PCS (x axis) and predicted PCS (y axis), labels both axes in ppm, uses blue circles and a red identity line, and sets both axis limits from the combined experimental and predicted values. The function prints the susceptibility tensor and magnetic multipole-centre location. It does not report a numeric fit metric or save a result file.

The distinctive model choice is the explicit [0, 1, 2] selector passed to `ilpcs`; unlike `tm_1igv_point_fit.m`, this call is not described as a point-electron fit.
