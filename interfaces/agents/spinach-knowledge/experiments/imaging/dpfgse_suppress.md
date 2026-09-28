# experiments/imaging/dpfgse_suppress.m

- Signature: `fid=dpfgse_suppress(spin_system,parameters,H,R,K,G,F)`

## Purpose

DPFGSE signal suppression, based on Equation 3 from the paper by Stott et al. (https://doi.org/10.1006/jmre.1997.1110). Syntax: fid=dpfgse_suppress(spin_system,parameters,H,R,K,G,F)

## Physical / mathematical content

The sequence starts with a hard proton 90-degree pulse, then applies two selective soft pulses interleaved with hard 180-degree pulses and pairs of gradient periods. Evolution uses `L=H+F+1i*R+1i*K`, with gradient periods adding `parameters.g_amp(1)*G{1}` or `parameters.g_amp(2)*G{1}`.

## Numerical / algorithmic content

Hard pulses and gradient periods propagate the state with `step`; selective pulses use `shaped_pulse_af` with the `'expv'` method and frequencies shifted by `parameters.offset`. The FID is computed by `evolution` from `parameters.coil` at intervals of `1/parameters.sweep`.

## Parameters / inputs

- parameters.g_amp -amplitudes of the two gradients, T/m
- parameters.g_dur -gradient duration, seconds
- parameters.rf_frq_list -soft pulse parameters that will
- parameters.rf_amp_list be passed to shaped_pulse_af
- parameters.rf_dur_list function
- parameters.rf_phi
- parameters.max_rank
- parameters.sweep -detection sweep width, Hz
- parameters.npoints -number of points in the fid
- parameters.offset -transmitter and receiver offset, Hz

## Outputs

- fid -free induction decay of what is effectively a
- 1D pulse-acquire NMR experiment
- Notes: at least a hundred points are required in the spatial
- dimension; increase until the answer stops changing.

## Implementation structure

The function validates required fields with the local `grumble` helper, constructs proton `Lx` and `Ly` operators across the spatial grid, and applies the pulse–gradient sequence to `parameters.rho0`. It then returns the coil-observed FID.
