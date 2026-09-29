# examples/nmr_overtone/cpmas_glycine_match_1.m

- MATLAB implementation: [examples/nmr_overtone/cpmas_glycine_match_1.m](https://github.com/IlyaKuprov/Spinach/blob/main/examples/nmr_overtone/cpmas_glycine_match_1.m)

- Signature: `cpmas_glycine_match_1()`

## What it models

This is a Spinach simulation wrapper for cross-polarisation from `1H` to a `14N` overtone transition in glycine under magic-angle spinning (MAS). That description comes from the source header. The calculation is not a fit to or processing of an acquired spectrum: the file contains no experimental-data read, and it builds its spin system from inline parameters before calling `singlerot(spin_system,@overtone_cp,parameters,'qnmr')`.

The source comment attributes the glycine quadrupolar tensor to O'Dell and Ratcliffe, [DOI 10.1016/j.cplett.2011.08.030](https://doi.org/10.1016/j.cplett.2011.08.030). The source converts `C_q=1.18e6` Hz (1.18 MHz), asymmetry `eta_q=0.53`, and spin `I=1` into the nitrogen quadrupolar interaction with `eeqq2nqi`. These are literature-motivated model inputs, not a measured spectrum loaded by this example.

## Spin system and approximations

The isotope labels are `14N` and `1H`; the field input is `sys.magnet=14.10220742`. The source also sets `inter.zeeman.scalar={32.4,0.0}` and coordinates `[0,0,0]` and `[0,0,1.00]`; it does not annotate units for these literals in this wrapper. No separate N–H tensor is assigned by hand: `create` generates the point-dipolar coupling from the 1 Å coordinates. A second dipolar tensor would double-count it. The RF/sequence calculation is delegated to `@overtone_cp`; its internal pulse program is not present in this file. No gradient list or laboratory pulse-acquire schedule is defined in this wrapper.

Relaxation is configured as `{'damp'}`, with `rlx_keep='diagonal'`, zero equilibrium, and `damp_rate=300` (the wrapper does not state a unit for that rate). The basis is spherical-tensor Liouville space with no approximation, and the spectrum calculation uses `max_rank=7`. The rotor-axis vector is `[sqrt(2/3),0,sqrt(1/3)]`; the RF/operator-state mixing angle is set to `atan(sqrt(2))`, the magic angle.

## RF sweep and simulated readout

The wrapper fixes the encoded rotor-rate parameter at `-19840` and selects the rough powder grid `rep_2ang_200pts_oct`. It sets a `44e3-52e3` frequency sweep (44-52 kHz on the displayed axis), 256 points, and 256-point zero filling. The initial state is an oriented `1H` state and the receive operator an oriented `14N` state.

The `14N` RF frequency is `48e3` (48 kHz). The contact-duration input is `1e-4` s. For each of 15 settings, the source scans the `1H` RF nutation-frequency input from 25 to 39 kHz in equal steps and supplies the two-channel power vector as `2*pi*[55e3,rf_powers(n)]/sin(theta)`. Each call produces one simulated spectrum; the code plots its real part in a 1-by-15 panel layout and fixes the display window to x=44-52 and y=-9e-4 to 9e-4. The y limits are plotting choices, not measured signal bounds. These axes and traces are calculation outputs, not reported experimental measurements. The header's “minutes” is a runtime estimate, not a timing measured in this task.
