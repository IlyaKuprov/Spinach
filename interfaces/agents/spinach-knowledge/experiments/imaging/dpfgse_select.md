# experiments/imaging/dpfgse_select.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/experiments/imaging/dpfgse_select.m`
- Signature: `fid=dpfgse_select(spin_system,parameters,H,R,K,G,F)`
- Total lines: 186

## Purpose

DPFGSE signal selection, based on Equation 3 from the paper by Stott et al. (https://doi.org/10.1006/jmre.1997.1110). Syntax: fid=dpfgse_select(spin_system,parameters,H,R,K,G,F)

## Physical / mathematical content

- Imaging sequence implementations. They build spatially resolved Liouvillians that include gradients, slice-selection RF terms, diffusion, and acquisition operators.

## Numerical / algorithmic content

- Time propagation is explicit. In Spinach this usually means repeated application of matrix exponentials or propagator factorizations to density operators or state vectors in Hilbert/Liouville/Fokker-Planck space.
- The file contains an explicit `grumble(...)` validator, which is Spinach convention for front-loading dimension, type, and regime checks before expensive linear-algebra work begins.
- The file also defines local helper function(s): `grumble()`. This usually means the public entry point is supported by tightly coupled validation or helper logic kept private to the file.

## Code-derived implementation details

### Comment-guided execution stages

- Lines 40-41: Check consistency; implemented by `grumble(parameters)`.
- Lines 43-44: Compose Liouvillian; implemented by `L=H+F+1i*R+1i*K`.
- Lines 46-47: Get pulse operators; implemented by `Lp=operator(spin_system,'L+','1H')`.
- Lines 51-52: Hard 90 on everything; implemented by `rho=step(spin_system,Ly,parameters.rho0,pi/2)`.
- Lines 54-55: Gradient; implemented by `rho=step(spin_system,L+parameters.g_amp(1)*G{1},rho,parameters.g_dur)`.
- Lines 57-60: Soft 180 on user-specified frequency; implemented by `rho=shaped_pulse_af(spin_system,L,Lx,Ly,rho,parameters.rf_frq_list-parameters.offset, parameters.rf_amp_list,parameters.rf_dur_list, parameters.rf_phi,parameters.max_ran…`.
- Lines 62-63: Gradients; implemented by `rho=step(spin_system,L+parameters.g_amp(1)*G{1},rho,parameters.g_dur)`.
- Lines 71-72: Gradient; implemented by `rho=step(spin_system,L+parameters.g_amp(2)*G{1},rho,parameters.g_dur)`.
- Lines 74-76: Run the evolution and watch the coil state; implemented by `fid=evolution(spin_system,L,parameters.coil,rho, 1/parameters.sweep,parameters.npoints-1,'observable')`.

### Key state/data transformations

- Lines 44: computes `L` using `L=H+F+1i*R+1i*K`.
- Lines 47: computes `Lp` using `Lp=operator(spin_system,'L+','1H')`.
- Lines 48: computes `Lx` using `Lx=kron(speye(prod(parameters.npts)),(Lp+Lp')/2)`.
- Lines 49: computes `Ly` using `Ly=kron(speye(prod(parameters.npts)),(Lp-Lp')/2i)`.
- Lines 52: computes `rho` using `rho=step(spin_system,Ly,parameters.rho0,pi/2)`.
- Lines 75-76: computes `fid` using `fid=evolution(spin_system,L,parameters.coil,rho, 1/parameters.sweep,parameters.npoints-1,'observable')`.

### Local helper functions

- Line 81: `grumble()` — `function grumble(parameters)`.
  - Representative operation: `if ~isfield(parameters,'rho0')`.
  - Representative operation: `error('parameters.rho0 field must be present.')`.

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

- DPFGSE signal selection, based on Equation 3 from the paper by
- Stott et al. (https://doi.org/10.1006/jmre.1997.1110). Syntax:
- fid=dpfgse_select(spin_system,parameters,H,R,K,G,F)
- parameters.g_amp -amplitudes of the two gradients, T/m
- parameters.g_dur -gradient duration, seconds
- parameters.rf_frq_list -soft pulse parameters that will
- parameters.rf_amp_list be passed to shaped_pulse_af
- parameters.rf_dur_list function
- parameters.rf_phi
- parameters.max_rank
- parameters.sweep -detection sweep width, Hz
- parameters.npoints -number of points in the fid

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `grumble()`, `operator()`, `speye()`, `step()`, `shaped_pulse_af()`, `evolution()`, `isfield()`, `isvector()`, `any()`, `isscalar()`.
