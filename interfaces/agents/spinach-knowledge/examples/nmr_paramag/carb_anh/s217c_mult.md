# examples/nmr_paramag/carb_anh/s217c_mult.m

- MATLAB implementation: [examples/nmr_paramag/carb_anh/s217c_mult.m](https://github.com/IlyaKuprov/Spinach/blob/main/examples/nmr_paramag/carb_anh/s217c_mult.m)

- Signature: `s217c_mult()`

## Purpose

Fits PCS for the S217C tag site of human carbonic anhydrase II (hCA-II) with the multipolar distributed-paramagnet model. It is distinct from both a single point-electron fit and the voxelised density reconstruction used by `s217c_kuprov`. The cited study describes Tm3+-DOTA-M8-tagged mutants and 1H-15N HSQC PCS referenced against diamagnetic Lu3+-DOTA-M8 samples; it reports 364 unambiguous 1H and 15N assignments for S217C. The script loads `s217c_expt.mat` but does not identify which nuclei or assignments are included there.

## Fit and outputs

The call `ilpcs(xyz,expt_pcs,[0 1 2],[-23 -16 20])` fits multipole ranks 0, 1, and 2 from measured PCS at nuclear coordinates; the solver documentation says rank 0 is fixed by normalisation, rather than fitted. The initial paramagnetic-centre coordinate is [-23, -16, 20] Angstrom. The fitted output includes centre coordinate `mxyz`, susceptibility tensor `chi`, multipole coefficients (not retained by this script), and predicted PCS. Coordinates are in Angstrom, PCS in ppm, and the tensor in Angstrom^3. The script plots experimental against predicted PCS and displays the tensor and centre location.

## Scope and limits

No multipole coefficients, actual fitted values, or residual statistics are saved. The script gives no magnetic field or acquisition temperature; the cited article narrative describes room-temperature solution measurements but no numeric field or temperature. Its MAT input does not expose the exact nuclear subset in the MATLAB source.

## Sources

- [MATLAB example source](https://github.com/IlyaKuprov/Spinach/blob/main/examples/nmr_paramag/carb_anh/s217c_mult.m)
- [Multipolar PCS fitter and units](https://github.com/IlyaKuprov/Spinach/blob/main/experiments/pseudocon/ilpcs.m)
- [Suturina et al., Chemical Science 8, 2751-2757 (2017), DOI: 10.1039/c6sc03736d](https://doi.org/10.1039/c6sc03736d)
- [PCS analysis tutorial](https://spindynamics.org/wiki/index.php?title=Pseudocontact_shift_analysis)
