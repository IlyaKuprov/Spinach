# examples/nmr_overtone/cpmas_valine_accum.m

- MATLAB implementation: [examples/nmr_overtone/cpmas_valine_accum.m](https://github.com/IlyaKuprov/Spinach/blob/main/examples/nmr_overtone/cpmas_valine_accum.m)

- Signature: `cpmas_valine_accum()`

## What it models

The source header describes cross-polarisation from protons to the `14N` overtone transition in N-acetylvaline under MAS, using Fokker-Planck formalism. The wrapper makes an accumulation profile by changing the RF contact duration and simulating a spectrum at each setting. The calculation is not based on a loaded measured spectrum: all model inputs are specified in the MATLAB file, which calls `singlerot(spin_system,@overtone_cp,parameters,'qnmr')`. Details of the pulse program inside `@overtone_cp` are not in this wrapper and are not inferred here. No gradient list or laboratory pulse-acquire schedule is defined here.

The source comment attributes the valine quadrupolar tensor data to [DOI 10.1039/c4cp03994g](https://doi.org/10.1039/c4cp03994g). The code converts `C_q=3.21e6` Hz (3.21 MHz), asymmetry `eta_q=0.27`, and `I=1` using `eeqq2nqi`. This is provenance for a numerical model input, not an external experimental file read by this script.

## Spin system and relaxation

The isotope labels are `14N` and `1H`, and the field input is `sys.magnet=14.10220742`. The source supplies Zeeman eigenvalue arrays `[57.5,81.0,227.0]` and `[0,0,0]`, with Euler angles `[-90,-90,-17]` degrees (the code multiplies by `pi/180`); the wrapper does not label the eigenvalue units. Coordinates are `[0,0,0]` and `[1,0,0]`, without a unit stated here. Damping relaxation retains diagonal terms, uses zero equilibrium, and sets `damp_rate=2000` (unit not stated in the wrapper).

The basis is spherical-tensor Liouville space without approximation. The source header identifies the approach as Fokker-Planck, while the visible wrapper sets `bas.formalism='sphten-liouv'` and disables `krylov` and `trajlevel`; no additional pulse-program details are specified. It sets `max_rank=9`, rotor-rate input `-19840`, and powder grid `rep_2ang_800pts_sph`.

## Contact-time sweep and readout

The frequency axis spans 70-105 kHz, with 256 points and 256-point zero filling. The initial state is for `1H`; receiver and overtone operators are for `14N`. The script uses `rf_frq=86.30e3` (86.30 kHz) and an RF power vector `2*pi*[55.0e3,35.1e3]/sin(theta)`, where `theta=atan(sqrt(2))`. In a 10-panel loop it sets `rf_dur=1e-5*n` s, so the encoded contact-duration values are 10, 20, ..., 100 microseconds.

Each panel plots the real part of a simulated spectrum. The displayed intensity limits, from `-1.7569e-3` to `1.7569e-3`, are plotting limits set by the source, not measured bounds. The source header's “hours” is an estimate; this task did not run the example or observe experimental data.
