# experiments/imaging/spin_echo.m

[MATLAB source](https://github.com/IlyaKuprov/Spinach/blob/main/experiments/imaging/spin_echo.m) · [Spinach Wiki](https://spindynamics.org/wiki/index.php?title=spin_echo.m)

## Purpose

`spin_echo` is an imaging-context spin-echo acquisition. It applies a hard 90° pulse, evolves under the supplied Liouvillian and the first spatial-gradient operator, applies a hard 180° pulse, then detects during a second gradient-evolution interval. This is a parameterised sequence, not a report of a measured experiment or a run-verified image.

## Inputs and units

The function is called by `imaging()` with `H`, `R`, `K`, `G`, and `F`. It assembles `L=H+F+1i*R+1i*K`. Its required sequence fields are:

- `parameters.g_amp`: gradient amplitude in T/m.
- `parameters.g_step_dur`: propagation-step duration in seconds.
- `parameters.g_n_steps`: positive integer number of steps in the first gradient interval; detection uses twice this count.
- `parameters.spins`: nonempty cell array of spin-name strings; the first entry selects the nucleus for the `L+` pulse operator.
- `parameters.rho0`, `parameters.coil`, and `parameters.npts`: initial state, detection observable, and spatial discretisation supplied by the imaging setup.

`G{1}` is used for both gradient intervals. The pulse operator is the y component formed from the selected `L+` operator and replicated over `npts`; the sequence does not select coherence orders or project a user-specified detection state.

## Detection and return

The return value `fid` is the observable detected with `parameters.coil` over `2*parameters.g_n_steps` propagation steps, each of duration `g_step_dur`, under `L+g_amp*G{1}`. It is the time-domain echo signal; the imaging grid is used internally by the spatial-gradient evolution. The function returns the signal only, without explicit time or spatial coordinate vectors.
