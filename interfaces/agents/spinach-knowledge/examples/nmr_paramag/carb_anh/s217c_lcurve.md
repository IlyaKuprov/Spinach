# examples/nmr_paramag/carb_anh/s217c_lcurve.m

- MATLAB implementation: [examples/nmr_paramag/carb_anh/s217c_lcurve.m](https://github.com/IlyaKuprov/Spinach/blob/main/examples/nmr_paramag/carb_anh/s217c_lcurve.m)

- Signature: `s217c_lcurve()`

## Purpose

Selects a Tikhonov regularisation parameter for a distributed PCS reconstruction of the S217C tag site of human carbonic anhydrase II (hCA-II). This is a regularisation diagnostic, not a point fit and not the final reconstruction. The cited study describes Tm3+-DOTA-M8-tagged mutants and 1H-15N HSQC PCS referenced against diamagnetic Lu3+-DOTA-M8 samples; it reports 364 unambiguous 1H and 15N assignments for S217C. The script loads the dataset and effective susceptibility tensor from MAT files, without recording which assignment subset is present.

## L-curve calculation

The inverse-problem setup uses the Kuprov equation, box centre [-21.8, -18.4, 20.2] Angstrom, 50 Angstrom sides, 50 Angstrom margins, confinement radii [2, 12] Angstrom, zero sharpening weight, and GPU execution. For each of 30 values from 10^-2 through 10^2, it calls `ipcs(parameters,64,lambda)` in a parallel loop, records least-squares error and regularisation output, and divides the latter by lambda. It calls `lcurve(lam,err,reg,'log')` and prints the suggested smoothing parameter. It does not save a density reconstruction or fitted tensor.

The solver documents coordinates and confinement distances in Angstrom, PCS in ppm, and susceptibility in Angstrom^3. The MATLAB source contains no exact magnetic field or acquisition temperature; the cited article narrative describes room-temperature solution measurements without a numeric temperature or field. The paper reports 1H and 15N shifts, but the loaded MAT file's exact nuclear subset is not specified by this script.

## Sources

- [MATLAB example source](https://github.com/IlyaKuprov/Spinach/blob/main/examples/nmr_paramag/carb_anh/s217c_lcurve.m)
- [Distributed PCS solver and units](https://github.com/IlyaKuprov/Spinach/blob/main/experiments/pseudocon/ipcs.m)
- [L-curve routine](https://github.com/IlyaKuprov/Spinach/blob/main/kernel/utilities/lcurve.m)
- [Suturina et al., Chemical Science 8, 2751-2757 (2017), DOI: 10.1039/c6sc03736d](https://doi.org/10.1039/c6sc03736d)
- [PCS analysis tutorial](https://spindynamics.org/wiki/index.php?title=Pseudocontact_shift_analysis)
