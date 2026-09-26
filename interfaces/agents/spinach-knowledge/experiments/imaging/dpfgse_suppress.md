# experiments/imaging/dpfgse_suppress.m

- Signature: `fid=dpfgse_suppress(spin_system,parameters,H,R,K,G,F)`

## Purpose

DPFGSE signal suppression, based on Equation 3 from the paper by Stott et al. (https://doi.org/10.1006/jmre.1997.1110). Syntax: fid=dpfgse_suppress(spin_system,parameters,H,R,K,G,F)

## Physical / mathematical content

- Imaging sequence implementations. They build spatially resolved Liouvillians that include gradients, slice-selection RF terms, diffusion, and acquisition operators.

## Numerical / algorithmic content

- Time propagation is explicit. In Spinach this usually means repeated application of matrix exponentials or propagator factorizations to density operators or state vectors in Hilbert/Liouville/Fokker-Planck space.

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

- DPFGSE signal suppression, based on Equation 3 from the paper by
- Stott et al. (https://doi.org/10.1006/jmre.1997.1110). Syntax:
- fid=dpfgse_suppress(spin_system,parameters,H,R,K,G,F)
- parameters.g_amp -amplitudes of the two gradients, T/m
- parameters.g_dur -gradient duration, seconds
- parameters.rf_frq_list -soft pulse parameters that will
- parameters.rf_amp_list be passed to shaped_pulse_af
- parameters.rf_dur_list function
- parameters.rf_phi
- parameters.max_rank
- parameters.sweep -detection sweep width, Hz
- parameters.npoints -number of points in the fid
