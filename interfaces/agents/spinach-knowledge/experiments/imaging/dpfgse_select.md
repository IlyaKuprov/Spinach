# experiments/imaging/dpfgse_select.m

- Signature: `fid=dpfgse_select(spin_system,parameters,H,R,K,G,F)`

## Purpose

DPFGSE signal selection, based on Equation 3 from the paper by Stott et al. (https://doi.org/10.1006/jmre.1997.1110). Syntax: fid=dpfgse_select(spin_system,parameters,H,R,K,G,F)

## Physical / mathematical content

The sequence applies a hard 90° pulse, then two soft 180° pulses, each bracketed by a pair of equal-amplitude gradient periods. Evolution uses `L=H+F+1i*R+1i*K`, with each gradient period adding `parameters.g_amp(i)*G{1}`.

## Numerical / algorithmic content

The hard pulse and gradient periods use `step`; each soft pulse uses `shaped_pulse_af` with RF frequencies shifted by `parameters.offset` and the `'expv'` method. The FID is acquired with `evolution` in `'observable'` mode at intervals of `1/parameters.sweep`.

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

The function validates the parameter fields with the local `grumble` function, constructs spatially extended `Lx` and `Ly` pulse operators, and propagates `parameters.rho0` through the pulse–gradient sequence. It then observes the resulting state through `parameters.coil` to return `fid`.
