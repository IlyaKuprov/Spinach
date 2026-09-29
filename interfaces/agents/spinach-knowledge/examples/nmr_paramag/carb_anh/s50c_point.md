# examples/nmr_paramag/carb_anh/s50c_point.m

- Signature: `s50c_point()`
- Source: [s50c_point.m](https://github.com/IlyaKuprov/Spinach/blob/main/examples/nmr_paramag/carb_anh/s50c_point.m)

## Purpose

Fit a point paramagnetic centre to experimental pseudocontact shifts (PCS) for the S50C mutant dataset of human carbonic anhydrase II. The source cites the [study](https://doi.org/10.1039/c6sc03736d) and a [PCS tutorial](https://spindynamics.org/wiki/index.php?title=Pseudocontact_shift_analysis).

## Model and inputs

- This is a point-centre inverse PCS fit, not a distributed electron-density calculation.
- Loads `expt_pcs` and nuclear coordinates `xyz` from `s50c_expt.mat`; the source does not state coordinate units or enumerate the nuclei/isotopes in this MAT-file.
- Calls `ippcs(xyz,[-27.0 13.0 18.0],expt_pcs)`. This vector is passed as the second argument; the example does not label its physical meaning or units.
- The example gives no magnetic-field or temperature value.

## Outputs and limits

The call returns `mxyz`, `chi`, and `pred_pcs`: the fitted point-electron location, susceptibility tensor, and predicted PCS, respectively. It plots experimental against predicted PCS with axes labelled in ppm and displays `chi` and `mxyz`. The source provides no numerical fit result in advance and does not state units for the tensor or returned coordinates. It does not simulate a spectrum.
