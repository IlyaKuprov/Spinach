# kernel/optimcon/tgrape.m

- Signature: `[fidelity,grad]=tgrape(spin_system,drift,controls,waveform,...`

`tgrape` evaluates a GRAPE control-sequence fidelity and its gradient with respect to slice durations. It supports only `sphten-liouv` and `zeeman-liouv` formalisms.

`drift` is a square Liouvillian; `controls` is a cell array of same-size square operators. `rho_init` and `rho_targ` are compatible column vectors. `waveform` is a real array with one row per control and one column per time slice; coefficients are in rad/s. `dt_grid` is a finite real column vector with one duration per slice. The positive scalar `time_unit`, in seconds, scales those durations.

`fidelity` is the real part of the target-state overlap. `grad` gives its derivative with respect to each `dt_grid` entry.

Source: https://spindynamics.org/wiki/index.php?title=tgrape.m