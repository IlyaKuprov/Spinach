# experiments/imaging/dpfgse_select.m

- Signature: `fid=dpfgse_select(spin_system,parameters,H,R,K,G,F)`
- Canonical MATLAB source: [`experiments/imaging/dpfgse_select.m`](https://github.com/IlyaKuprov/Spinach/blob/main/experiments/imaging/dpfgse_select.m)

## Contract and physical sequence

This is a simulated one-dimensional DPFGSE signal-selection / pulse-acquire experiment, not an MRI image or a DNP calculation. The source attributes the sequence to Equation 3 of Stott et al. It forms `L=H+F+1i*R+1i*K`, uses the first spatial gradient operator `G{1}`, and constructs proton `Lx` and `Ly` pulse operators expanded over `prod(parameters.npts)` spatial points.

The initial state `parameters.rho0` receives a hard proton pi/2 about `Ly`. The ordered train is: a `g_amp(1)` gradient; one shaped soft pi pulse; gradients `g_amp(1)` then `g_amp(2)`; a second shaped soft pi pulse; and a final `g_amp(2)` gradient. Each gradient period lasts `parameters.g_dur`. Both shaped pulses are applied by `shaped_pulse_af` at the listed frequencies shifted by `parameters.offset`; amplitude, duration, phase, and rank controls are passed through the corresponding parameter fields. This is the source's signal-selection variant; it reports no measured selection profile or signal.

## Parameters and units

- `parameters.g_amp`: two gradient amplitudes, T/m; the gradient durations are `parameters.g_dur`, seconds.
- `parameters.rf_frq_list`, `parameters.rf_amp_list`, `parameters.rf_dur_list`, `parameters.rf_phi`, and `parameters.max_rank`: shaped-pulse inputs forwarded to `shaped_pulse_af`; the source checks equal element counts for the frequency, amplitude, and duration lists and positive pulse durations.
- `parameters.rho0` and numeric `parameters.coil`: initial state and detection state.
- `parameters.npts`: positive-integer spatial-grid vector. The source notes that at least 100 spatial points are required and recommends increasing the grid until the answer stops changing; this is a resolution guideline, not a reported simulation result.
- `parameters.offset` and `parameters.sweep` are in Hz; `parameters.npoints` is the FID point count. Acquisition passes `1/sweep` seconds and `npoints-1` evolution steps to `evolution`, observing `parameters.coil`.

## References

- Stott et al., DPFGSE method, cited by the source: <https://doi.org/10.1006/jmre.1997.1110>
- Spinach Wiki: <https://spindynamics.org/wiki/index.php?title=dpfgse_select.m>
- Source: [`experiments/imaging/dpfgse_select.m`](https://github.com/IlyaKuprov/Spinach/blob/main/experiments/imaging/dpfgse_select.m)
