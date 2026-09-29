# examples/nmr_paramag/carb_anh/s217c_kuprov.m

- MATLAB implementation: [examples/nmr_paramag/carb_anh/s217c_kuprov.m](https://github.com/IlyaKuprov/Spinach/blob/main/examples/nmr_paramag/carb_anh/s217c_kuprov.m)

- Signature: `s217c_kuprov()`

## Purpose

Reconstructs a spatially distributed PCS source for the S217C tag site of human carbonic anhydrase II (hCA-II), rather than fitting one electron point. The cited study describes Tm3+-DOTA-M8-tagged mutants and 1H-15N HSQC PCS referenced against diamagnetic Lu3+-DOTA-M8 samples; it reports 364 unambiguous 1H and 15N assignments for S217C. The script loads `s217c_expt.mat`; that source does not say which assignment subset the MAT file contains.

## Density fit and outputs

The code loads a precomputed effective susceptibility tensor from `s217c_chi_eff.mat` and configures `ipcs` with the Kuprov equation. The reconstruction box is centred at [-21.8, -18.4, 20.2] Angstrom and has 25 Angstrom sides; the confinement radii are 3 and 12 Angstrom, sharpening weight is 1, and six margins are 50 Angstrom each. PCS, nuclear positions, all-atom coordinates, and the tensor are passed to the solver; `gpu=true()` enables its GPU option.

It solves successively on 64^3, 128^3, and 256^3 grids with regularisation parameter 0.50, feeding each source cube into the next grid as the guess. It then recalculates an effective susceptibility tensor with `chi_eff(source_cube,ranges,xyz,expt_pcs)` and displays that tensor. The source does not report a density centroid, tensor values, or PCS fit statistics; those are data- and run-dependent.

## Scope and limits

The `ipcs` interface documents nuclear and atom coordinates, box coordinates, and confinement distances in Angstrom, PCS in ppm, and susceptibility in Angstrom^3. Neither this script nor the cited article narrative specifies an exact magnetic field or numeric temperature; the article describes solution measurements at room temperature. No single electron coordinate is fitted by this density workflow.

## Sources

- [MATLAB example source](https://github.com/IlyaKuprov/Spinach/blob/main/examples/nmr_paramag/carb_anh/s217c_kuprov.m)
- [Distributed PCS solver and units](https://github.com/IlyaKuprov/Spinach/blob/main/experiments/pseudocon/ipcs.m)
- [Effective susceptibility fitter](https://github.com/IlyaKuprov/Spinach/blob/main/experiments/pseudocon/chi_eff.m)
- [Suturina et al., Chemical Science 8, 2751-2757 (2017), DOI: 10.1039/c6sc03736d](https://doi.org/10.1039/c6sc03736d)
- [PCS analysis tutorial](https://spindynamics.org/wiki/index.php?title=Pseudocontact_shift_analysis)
