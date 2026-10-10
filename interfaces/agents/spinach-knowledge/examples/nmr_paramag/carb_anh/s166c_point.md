# examples/nmr_paramag/carb_anh/s166c_point.m

- MATLAB implementation: [examples/nmr_paramag/carb_anh/s166c_point.m](https://github.com/IlyaKuprov/Spinach/blob/main/examples/nmr_paramag/carb_anh/s166c_point.m)

- Signature: `s166c_point()`

## Purpose

Fits pseudocontact shifts (PCS) for the S166C tag site of human carbonic anhydrase II (hCA-II) with a single point-electron model. The cited study describes Tm3+-DOTA-M8-tagged hCA-II mutants and PCS from 1H-15N HSQC spectra referenced against diamagnetic Lu3+-DOTA-M8 samples; it reports 397 unambiguous 1H and 15N assignments for S166C. The script loads the prepared arrays from `s166c_expt.mat`; it does not identify which subset of those assignments is in that file.

## Fit and outputs

The call `ippcs(xyz,[-15 -3 -10],expt_pcs)` uses the nuclear positions and measured shifts to fit a point location `mxyz`, susceptibility tensor `chi`, and predicted shifts. The starting electron coordinate is [-15, -3, -10] Angstrom; the solver documents nuclear coordinates and PCS in Angstrom and ppm, respectively. Its susceptibility tensor is in Angstrom^3. The script plots experimental against predicted PCS and displays the fitted tensor and point location; it does not write a result file.

## Scope and limits

The MATLAB source specifies no magnetic field or acquisition temperature. The paper describes room-temperature solution measurements but the cited narrative gives no numerical field or temperature. The script contains no observed fit results; they depend on the loaded MAT file. Coordinate and tensor values printed at runtime are outputs, not constants supplied here.

## Sources

- [MATLAB example source](https://github.com/IlyaKuprov/Spinach/blob/main/examples/nmr_paramag/carb_anh/s166c_point.m)
- [Point PCS fitter and units](https://github.com/IlyaKuprov/Spinach/blob/main/experiments/pseudocon/ippcs.m)
- [Suturina et al., Chemical Science 8, 2751-2757 (2017), DOI: 10.1039/c6sc03736d](https://doi.org/10.1039/c6sc03736d)
- [PCS analysis tutorial](https://spindynamics.org/wiki/index.php?title=Pseudocontact_shift_analysis)
