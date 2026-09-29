# experiments/esr_dipolar/deer_3p_hard_echo.m

This is an echo-window sampler for the three-pulse DEER experiment, not a routine that constructs a complete DEER trace. The source describes its diagnostic use: locating an echo that may be narrow in simulation and, for high-spin electrons, displaced from its expected position. It uses the supplied hard-pulse operators as-is; their spin or transition selectivity is determined by the caller.

## Sequence and detection

Starting from `parameters.rho0`, the function applies `parameters.ex_prob` at `pi/2`, evolves for `parameters.tb`, applies `parameters.ex_pump` at `pi`, evolves for `parameters.ta-parameters.tb`, then applies `parameters.ex_prob` at `pi`. After that pulse it evolves for `parameters.ta-parameters.tc/2` and records the observable with `parameters.coil` across an interval of `parameters.tc` using `parameters.nsteps+1` samples including the initial observation, spaced by `parameters.tc/parameters.nsteps`. Thus the third pulse is at time `parameters.ta`, and the sampled window is centred at `2*parameters.ta` relative to the sequence start. Here `parameters.ta` is the first-to-third-pulse time, `parameters.tb` the first-to-second-pulse time, and `parameters.tc` the acquisition-window width; all are in seconds.

The required fields are `parameters.ex_prob` (probe operator), `parameters.ex_pump` (pump operator), `parameters.rho0` (initial state), `parameters.coil` (detection state), `parameters.ta`, `parameters.tb`, `parameters.tc`, and `parameters.nsteps` (positive integer propagation-step count, yielding `nsteps+1` echo samples). The context supplies `H`, `R`, and `K` as same-sized matrices; they form `L=H+1i*R+1i*K`. The source checks the timing constraints `parameters.tb<=parameters.ta` and `parameters.tc/2<=parameters.ta`.

The return value `echo` is the sampled signal with the requested number of points. The source does not state its array orientation, numeric units, or whether a caller should interpret the signal as real or complex. It specifies this three-pulse sequence only; it does not define a CPMG/CP, Bruker, or four-pulse timing scheme.

Source: https://github.com/IlyaKuprov/Spinach/blob/main/experiments/esr_dipolar/deer_3p_hard_echo.m
