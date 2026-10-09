# examples/nmr_overtone/cpmas_glycine_match_2.m

- MATLAB implementation: [examples/nmr_overtone/cpmas_glycine_match_2.m](https://github.com/IlyaKuprov/Spinach/blob/main/examples/nmr_overtone/cpmas_glycine_match_2.m)

- Signature: `cpmas_glycine_match_2()`

## What it models

This wrapper simulates cross-polarisation from `1H` to the `14N` overtone transition in glycine under MAS, as described in the source header. It varies sample spinning rate and proton RF nutation frequency to map a Hartmann-Hahn matching profile over a rough powder grid. It does not load a measured spectrum: the spin system and acquisition parameters are defined in the MATLAB file and passed to `singlerot(spin_system,@overtone_cp,localpar,'qnmr')`.

The source attributes the glycine quadrupolar tensor to O'Dell and Ratcliffe, [DOI 10.1016/j.cplett.2011.08.030](https://doi.org/10.1016/j.cplett.2011.08.030). It encodes `C_q=1.18e6` Hz (1.18 MHz), `eta_q=0.53`, and `I=1` through `eeqq2nqi`; these literature-motivated values are model inputs, not a measured result from this program.

## Spin system and model settings

The isotopes are `14N` and `1H`, with `sys.magnet=14.10220742`. The wrapper also specifies `inter.zeeman.scalar={32.4,0.0}` and coordinates `[0,0,0]` and `[0,0,1.00]`; it does not label units for those literals. No distinct N–H tensor is assigned by hand: `create` generates the point-dipolar coupling from the 1 Å coordinates. A second dipolar tensor would double-count it. The source delegates pulse/sequence execution to `@overtone_cp` and does not expose that function's internal pulse details. No gradient list or laboratory pulse-acquire schedule is defined in this wrapper.

Relaxation is configured as damping with diagonal relaxation terms retained, zero equilibrium, and `damp_rate=1000` (unit not stated in this wrapper). The basis is spherical-tensor Liouville space without approximation; `max_rank=5`. The rotor rate is set per point, using negative values corresponding to the 20-90 kHz sweep. The rotor grid is `rep_2ang_200pts_oct`.

## Two-dimensional parameter sweep and output

The wrapper samples 50 proton RF settings from 10 to 200 kHz and 50 sample spinning rates from 20 to 90 kHz. For each pair it forms `rf_pwr=2*pi*[55e3,rf_powers(n)]/sin(theta)`; the rate is `-spin_rates(k)`, and the nitrogen RF-frequency input follows `8e3-2*localpar.rate`. It sets a 100-microsecond contact-duration input. Each simulated spectrum uses 256 points and 256-point zero filling, with a sweep extending 4 kHz on either side of the per-point RF frequency. The receive operator is for `14N`, and the initial state is for `1H`.

For each condition the program reduces the real spectrum to `sum(real(spectrum))`, then plots the resulting values as an image: proton RF nutation frequency in kHz horizontally and sample spinning rate in kHz vertically. The image is a computed matching-profile map; this source contains no measured intensity array or displayed numerical optimum. The header's “hours” describes expected calculation time, not a measured runtime here.
