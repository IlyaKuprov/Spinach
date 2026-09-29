# experiments/spen/psyche.m

Canonical source: https://github.com/IlyaKuprov/Spinach/blob/main/experiments/spen/psyche.m
Spinach Wiki: https://spindynamics.org/wiki/index.php?title=psyche.m

`fid=psyche(spin_system,parameters,H,R,K,G,F)` models the PSYCHE pure-shift NMR sequence and returns a two-dimensional FID. The imaging context supplies the Fokker–Planck operators; the source combines them as `L=H+F+1i*R+1i*K` and extends the spin pulse operators across the spatial grid.

## Coherence pathway, chirps, and acquisition

From `rho0`, a hard 90-degree pulse about `Lx` precedes the first half of F1 evolution. The sequence selects `+1` coherence, evolves for `delta` under `L+g_amp*G{1}`, applies a hard 180-degree pulse about `Lx`, and repeats that gradient-containing interval before selecting `-1` coherence. The PSYCHE element is a pair of chirps under the same gradient-containing generator: the first uses quadratures `{Cx,+Cy}`; after selection of zero coherence, the second uses `{Cx,-Cy}`. The pathway is then selected at `+1` and the second half of F1 is propagated in refocusing mode.

The waveform comes from `chirp_pulse(pulsenpoints,duration,bandwidth,smfactor,chirptype)`. The source normalises both components by `max(Cx)` and scales them by `2*pi*rfbeta`. For `chirptype='saltire'`, `rfbeta=(beta/360)*sqrt(2*bandwidth/duration)`; for other supported chirp types it computes `q_beta=-(2*log(cosd(beta)/2+1/2))/pi` and `rfbeta=sqrt(duration*bandwidth*q_beta/(2*pi))/duration`. Each chirp lasts `duration`, divided equally across the waveform points.

F1 is represented by `npoints(1)` points: each half uses `timestep1/2` and `npoints(1)-1` propagation steps, first in trajectory mode and then in refocusing mode. F2 is observed with `coil` under `L`, using `timestep2` and `npoints(2)-1` steps. Thus `fid` carries the requested F1 and F2 axes with point counts `npoints(1)` and `npoints(2)`; it is a simulated FID, not a transformed spectrum or measured result.

## Required inputs

The function requires `rho0`, `coil`, one-element `spins`, `g_amp` (T/m), `dims` (m), `npts`, `npoints` (two integer axis lengths, each greater than one), `diff` (m^2/s), `pulsenpoints`, `duration` (s), `bandwidth` (Hz), `smfactor`, `chirptype`, `beta` (degrees), `timestep1` and `timestep2` (s), and `delta` (s). The checks require positive `dims`, `npts`, `pulsenpoints`, `duration`, `bandwidth`, and time steps; non-negative `diff` and `delta`; and finite real scalar gradient amplitude, smoothing factor, and flip angle. Supported `chirptype` values are `wurst`, `wurst-adaptive`, `smoothed`, `smoothed-adaptive`, `saltire`, and `saltire-adaptive`. The source requires `sphten-liouv` formalism; matrix inputs `H`, `R`, `K`, and `F` must have matching dimensions, and `G` must be a three-element cell array of numeric gradient operators.
