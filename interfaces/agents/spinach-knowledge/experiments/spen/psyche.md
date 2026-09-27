# experiments/spen/psyche.m

- Signature: `fid=psyche(spin_system,parameters,H,R,K,G,F)`

## Purpose

PSYCHE pure-shift NMR pulse sequence that returns a two-dimensional free induction decay.

## Parameters / inputs

- `parameters.rho0`: initial state.
- `parameters.coil`: detection state.
- `parameters.spins`: nuclei on which the sequence runs; a one-element cell array.
- `parameters.g_amp`: gradient amplitude (T/m).
- `parameters.dims`: sample size (m).
- `parameters.npts`: number of spatial grid discretisation points.
- `parameters.npoints`: two-element point-count vector `[F1 F2]`.
- `parameters.diff`: diffusion constant (m^2/s).
- `parameters.delta`: gradient evolution delay on either side of the hard 180-degree pulse (s).
- `parameters.timestep1`, `parameters.timestep2`: F1 and F2 evolution time steps (s).
- `parameters.pulsenpoints`: number of chirp waveform discretisation points.
- `parameters.duration`: duration of each PSYCHE chirp pulse (s).
- `parameters.bandwidth`: chirp sweep bandwidth around zero frequency (Hz).
- `parameters.smfactor`: chirp smoothing parameter; see `chirp_pulse.m`.
- `parameters.chirptype`: chirp waveform type; see `chirp_pulse.m`.
- `parameters.beta`: PSYCHE-element flip angle (degrees).
- `H`, `R`, `K`: Fokker–Planck Hamiltonian, relaxation superoperator, and kinetics superoperator, respectively, from the imaging context.
- `G`: three Fokker–Planck gradient superoperators from the imaging context.
- `F`: Fokker–Planck diffusion and flow superoperator from the context.

## Physical / numerical content

The sequence forms `L=H+F+1i*R+1i*K` and constructs spatially extended `Lx` and `Ly` pulse operators. It generates a chirp waveform with `chirp_pulse`, calculates its RF amplitude from `beta` using a separate formula for `chirptype='saltire'`, normalises the waveform, and scales both quadratures by `2*pi*rfbeta`.

A hard 90-degree pulse precedes the first half of F1 evolution. The sequence selects `+1` coherence, evolves for `delta` under `L+g_amp*G{1}`, applies a hard 180-degree pulse, and repeats that gradient evolution. It then selects `-1` coherence, applies the first gradient-assisted chirp with quadratures `{Cx,+Cy}`, selects `0` coherence, applies the second chirp with `{Cx,-Cy}`, and selects `+1` coherence. Each chirp is propagated with `shaped_pulse_xy` using `expv-pwc`. Refocused second-half F1 evolution is followed by F2 observable evolution using `parameters.coil`.

## Output

- `fid`: PSYCHE free induction decay as a 2D array.

## Consistency checks

Requires `sphten-liouv` formalism; matching-dimension numeric matrices `H`, `R`, `K`, and `F`; and a three-element cell array of numeric gradient operators `G`. Initial and detection states must be numeric column vectors matching `H`. `npoints` contains two integers greater than one; `npts` and `pulsenpoints` are positive integers. `dims`, `duration`, `bandwidth`, `timestep1`, and `timestep2` are positive finite real scalars; `diff` and `delta` are non-negative finite real scalars; `g_amp`, `smfactor`, and `beta` are finite real scalars. Supported `chirptype` values are `wurst`, `wurst-adaptive`, `smoothed`, `smoothed-adaptive`, `saltire`, and `saltire-adaptive`.

## Source attribution

- mohammadali.foroozandeh@chem.ox.ac.uk
- mariagrazia.concilio@sjtu.edu.cn
- <https://spindynamics.org/wiki/index.php?title=psyche.m>