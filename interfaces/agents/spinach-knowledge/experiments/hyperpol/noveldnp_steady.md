# experiments/hyperpol/noveldnp_steady.m

- Signature: `dnp=noveldnp_steady(spin_system,parameters,H,R,K)`
- Canonical MATLAB source: https://github.com/IlyaKuprov/Spinach/blob/main/experiments/hyperpol/noveldnp_steady.m

## Purpose and sequence

Computes a steady-state detected observable over microwave resonance offsets for the pulsed solid-effect or NOVEL DNP sequence. The input relaxation superoperator `R` is documented as thermalised to finite temperature. The function forms `L=H+1i*R+1i*K`; for each offset it adds the electron (L_z) offset term, applies the selected irradiation/contact block, then the unirradiated shot-spacing block, and solves the repeated-cycle steady state with a Newton solver. The output is the coil-state overlap with that steady state, not a measured polarisation.

For `flippulse=1`, the source applies an electron (X)-axis 90-degree pulse followed by a (-Y) microwave contact period. When `flipback=1`, it then applies a (-X) microwave flipback pulse of duration `pulse_dur`. For `flippulse=0`, it uses the contact period without the preparation pulse. `irr_powers`, `el_offs`, and `addshift` enter the generator multiplied by `2*pi`; use Hz-valued frequencies consistent with the source's stated Hz convention for microwave amplitude. Durations are in seconds.

## Inputs and output

`H`, `R`, and `K` are context-supplied matrices; `R` must represent finite-temperature relaxation. Required fields are `irr_powers` (non-negative microwave amplitude, Hz), `coil` (detection-state column vector), `contact_dur` (seconds), `shot_spacing` (seconds), `flippulse` (0 or 1), `flipback` (0 or 1), `addshift` (real scalar frequency shift), and `el_offs` (real offset array). If `flippulse=1`, a positive `pulse_dur` in seconds is required. `dnp` has the same shape as `parameters.el_offs`, with one complex-capable detected value per offset.



No gradients, spatial encoding, k-space, or FID are produced. The code-level 0/1 switch and 0/90/270-degree pulse axes are the available numerical examples; no parameter set or computed experimental result is provided in the source or baseline page.

## References

- NOVEL/solid-effect references retained from the source: https://doi.org/10.1016/0022-2364(88)90190-4 and https://doi.org/10.1063/1.5000528
- Spin Dynamics Wiki: https://spindynamics.org/wiki/index.php?title=noveldnp_steady.m
