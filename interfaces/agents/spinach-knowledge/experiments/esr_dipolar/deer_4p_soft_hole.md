# experiments/esr_dipolar/deer_4p_soft_hole.m

- MATLAB implementation: [experiments/esr_dipolar/deer_4p_soft_hole.m](https://github.com/IlyaKuprov/Spinach/blob/main/experiments/esr_dipolar/deer_4p_soft_hole.m)

## Purpose

This function diagnoses the frequency-selective effect of the soft pulses used in four-pulse DEER/PELDOR. The source describes it as a hypothetical experiment: each selected soft pulse is followed by an ideal hard π/2 pulse on all spins and time-domain detection. It is a pulse-diagnostic callback, not a DEER echo calculation.

## Interface and parameters

`fids=deer_4p_soft_hole(spin_system,parameters,H,R,K)`

The required parameter fields are `pulse_frq`, `pulse_pwr`, `pulse_dur`, `pulse_phi`, and `pulse_rnk` (four pulse values; respectively Hz, rad/s, seconds, radians, and Fokker–Planck ranks); `offset` (Hz), `sweep` (time-domain sweep width in Hz), `npoints` (FID points), `spins` (irradiated spin labels, normally `{'E'}`), `rho0` (initial state), `coil` (detection state), and `method` (`'expm'`, `'expv'`, or `'evolution'`). `H`, `R`, and `K` are same-sized matrices supplied by the context function. The implementation converts to Liouville representation when needed and requires `sphten-liouv` or `zeeman-liouv` formalism.

## Pulse and acquisition meaning

The routine evaluates four separate soft-pulse cases from the same `rho0`, using the four pulse parameter sets; these are not applied sequentially as one four-pulse train. It collects the four resulting states, applies the common hard π/2 step, and passes them to acquisition. The header describes that hard pulse as acting on all spins, while the implementation constructs its pulse operator from `spins{1}`; the supplied spin selection therefore matters. The declared output is four free-induction decays to apodise and Fourier-transform. The diagnostic description calls the time-domain detection infinite-bandwidth; the acquisition parameters still include `offset`, `sweep`, and `npoints`.

The output description does not state the FID array orientation or a sample-value convention. The source gives no pulse values or default acquisition settings. Powder averaging is done by the companion diagnostic driver, not by this function itself.

Source: https://spindynamics.org/wiki/index.php?title=deer_4p_soft_hole.m
