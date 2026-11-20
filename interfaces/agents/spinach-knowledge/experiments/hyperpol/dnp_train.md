# Repeated microwave blocks for pulsed DNP

Use `powder(spin_system,@dnp_train,parameters,'esr')` for static solid DNP in the electron rotating frame with nuclei in the laboratory frame. `crystal` evaluates an individual orientation. Use a Liouville basis and pass the initial state and receiver as matching column vectors.

Give `parameters.pulse_ops={Ex,Ey}`, two Cartesian RF amplitude rows in rad/s in `parameters.amplitudes`, a matching row of positive durations in seconds in `parameters.slice_durs`, and a positive integer `parameters.nloops`. For nutation frequencies in Hz and phases in degrees use `2*pi*freq.*cosd(phase)` and `2*pi*freq.*sind(phase)`. Drift, static relaxation, and kinetics remain active during every slice. The block repeats without a gap; include an explicit zero-RF slice for a delay.

One block propagator per orientation and detuning is built with the piecewise-constant shaped-pulse engine. It is reused for all repeats, so this approach is intended for systems small enough to store the Liouville propagator. The complex signal includes the pre-pulse point and one point per block. Average the complex signals before taking their real parts; the last point is the terminal transfer at the chosen offset.
