# experiments/imaging/dpfgse_suppress.m

- Signature: `fid=dpfgse_suppress(spin_system,parameters,H,R,K,G,F)`
- Canonical MATLAB source: [`experiments/imaging/dpfgse_suppress.m`](https://github.com/IlyaKuprov/Spinach/blob/main/experiments/imaging/dpfgse_suppress.m)

## Contract and physical sequence

This is a simulated one-dimensional DPFGSE signal-suppression / pulse-acquire experiment, not an MRI image or a DNP calculation. The source cites Equation 3 of Stott et al. It forms `L=H+F+1i*R+1i*K`, uses `G{1}` for spatial encoding, and constructs proton `Lx` and `Ly` pulse operators across `prod(parameters.npts)` spatial points.

The source sequence is: hard proton pi/2 about `Ly`; one `g_amp(1)` gradient; a soft shaped pi pulse at the listed frequencies shifted by `offset`; a hard pi rotation about `Ly` (the code passes angle `-pi`); gradients `g_amp(1)` then `g_amp(2)`; a second identical shaped soft pi pulse; another hard `-pi` about `Ly`; and a final `g_amp(2)` gradient. Every gradient lasts `g_dur`. It then observes the state through `parameters.coil`. This describes the suppression variant's implemented ordering; no measured suppression factor is claimed.

## Parameters and units

- `parameters.g_amp`: two gradient amplitudes in T/m; `parameters.g_dur` is a positive duration in seconds.
- `parameters.rf_frq_list`, `parameters.rf_amp_list`, `parameters.rf_dur_list`, `parameters.rf_phi`, and `parameters.max_rank`: shaped-pulse inputs passed to `shaped_pulse_af`; the source checks matching frequency/amplitude/duration list lengths and positive pulse durations.
- `parameters.rho0` and numeric `parameters.coil`: initial state and detected spin state.
- `parameters.npts`: positive-integer spatial-grid vector. The source recommends at least 100 points in the spatial dimension and increasing until the answer stops changing; this is a resolution guideline, not a result from a run.
- `parameters.offset` and `parameters.sweep` are in Hz; `parameters.npoints` is the FID point count. Acquisition passes `1/sweep` seconds and `npoints-1` evolution steps to `evolution` with the observable coil state.

## References

- Stott et al., DPFGSE method, cited by the source: <https://doi.org/10.1006/jmre.1997.1110>
- Spinach Wiki: <https://spindynamics.org/wiki/index.php?title=dpfgse_suppress.m>
- Source: [`experiments/imaging/dpfgse_suppress.m`](https://github.com/IlyaKuprov/Spinach/blob/main/experiments/imaging/dpfgse_suppress.m)
