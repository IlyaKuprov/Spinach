# experiments/nmr_solids/wise.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/experiments/nmr_solids/wise.m`
- Signature: `fid=wise(spin_system,parameters,H,R,K)`
- Total lines: 176

## Purpose

WISE (WIdeline SEparation) is a powder MAS heteronuclear correlation experiment. In the common 1H-13C implementation, molecular dynamics information is contained in 1H line shapes that are separated in the second dimension by 13C chemical shifts. Further information in:

## Physical / mathematical content

- Solid-state pulse sequence implementations. The core ingredients are anisotropic Hamiltonians, rotor synchronisation, cross-polarisation, recoupling/decoupling, and powder or rotor-stack propagation.
- Orientation or trajectory averaging is performed numerically, so grid design, weights, and integration error control matter directly to accuracy and runtime.

## Numerical / algorithmic content

- Time propagation is explicit. In Spinach this usually means repeated application of matrix exponentials or propagator factorizations to density operators or state vectors in Hilbert/Liouville/Fokker-Planck space.
- Numerical integration over angles or geometry is part of the implementation, so point placement and weights are as important as the local Hamiltonian calculations.
- The file contains an explicit `grumble(...)` validator, which is Spinach convention for front-loading dimension, type, and regime checks before expensive linear-algebra work begins.
- The file also defines local helper function(s): `grumble()`. This usually means the public entry point is supported by tightly coupled validation or helper logic kept private to the file.

## Code-derived implementation details

### Comment-guided execution stages

- Lines 50-51: Consistency enforcement; implemented by `grumble(spin_system,parameters,H,R,K)`.
- Lines 53-54: Wipe the state of 13C (pre-saturation); implemented by `[~,parameters.rho0]=decouple(spin_system,[],parameters.rho0,{'13C'})`.
- Lines 56-57: Build 1H and 13C control operators; implemented by `Hp=operator(spin_system,'L+',parameters.spins{1})`.
- Lines 63-64: Compose the Liouvillian; implemented by `L=H+1i*R+1i*K`.
- Lines 66-67: High-power 90-degree pulses on 1H along X (cos) and Y (sin); implemented by `L_hp_cos=L+2*pi*parameters.hi_pwr*Hx; L_hp_sin=L+2*pi*parameters.hi_pwr*Hy`.
- Lines 71-72: Get dwell times; implemented by `dw=1./parameters.sweep`.
- Lines 74-76: Run the F1 evolution; implemented by `rho_stack_cos=evolution(spin_system,L,[],rho_cos,dw(1), parameters.npoints(1)-1,'trajectory')`.
- Lines 80-81: CP contact time evolution generator (-Y on 1H, +X on 13C); implemented by `L_cp=L-2*pi*parameters.cp_pwr(1)*Hy+2*pi*parameters.cp_pwr(2)*Cx`.
- Lines 83-85: Run CP contact time evolution; implemented by `rho_stack_cos=evolution(spin_system,L_cp,[],rho_stack_cos, parameters.cp_dur,1,'final')`.
- Lines 89-90: Wipe and decouple protons for acquisition; implemented by `[L_dec,rho_stack_cos]=decouple(spin_system,L,rho_stack_cos,parameters.spins(1))`.
- Lines 93-95: Run the F2 evolution; implemented by `fid.cos=evolution(spin_system,L_dec,parameters.coil,rho_stack_cos, dw(2),parameters.npoints(2)-1,'observable')`.

### Key state/data transformations

- Lines 54: computes `[~,parameters.rho0]` using `[~,parameters.rho0]=decouple(spin_system,[],parameters.rho0,{'13C'})`.
- Lines 57: computes `Hp` using `Hp=operator(spin_system,'L+',parameters.spins{1})`.
- Lines 58: computes `Cp` using `Cp=operator(spin_system,'L+',parameters.spins{2})`.
- Lines 61: computes `Hx` using `Hx=(Hp+Hp')/2; Hy=(Hp-Hp')/2i; Cx=(Cp+Cp')/2`.
- Lines 64: computes `L` using `L=H+1i*R+1i*K`.
- Lines 67: computes `L_hp_cos` using `L_hp_cos=L+2*pi*parameters.hi_pwr*Hx; L_hp_sin=L+2*pi*parameters.hi_pwr*Hy`.
- Lines 68: computes `rho_cos` using `rho_cos=step(spin_system,L_hp_cos,parameters.rho0,1/(4*parameters.hi_pwr))`.
- Lines 69: computes `rho_sin` using `rho_sin=step(spin_system,L_hp_sin,parameters.rho0,1/(4*parameters.hi_pwr))`.
- Lines 72: computes `dw` using `dw=1./parameters.sweep`.
- Lines 75-76: computes `rho_stack_cos` using `rho_stack_cos=evolution(spin_system,L,[],rho_cos,dw(1), parameters.npoints(1)-1,'trajectory')`.
- Lines 77-78: computes `rho_stack_sin` using `rho_stack_sin=evolution(spin_system,L,[],rho_sin,dw(1), parameters.npoints(1)-1,'trajectory')`.
- Lines 81: computes `L_cp` using `L_cp=L-2*pi*parameters.cp_pwr(1)*Hy+2*pi*parameters.cp_pwr(2)*Cx`.
- Lines 90: computes `[L_dec,rho_stack_cos]` using `[L_dec,rho_stack_cos]=decouple(spin_system,L,rho_stack_cos,parameters.spins(1))`.
- Lines 91: computes `[~,rho_stack_sin]` using `[~,rho_stack_sin]=decouple(spin_system,[],rho_stack_sin,parameters.spins(1))`.
- Lines 94-95: computes `fid.cos` using `fid.cos=evolution(spin_system,L_dec,parameters.coil,rho_stack_cos, dw(2),parameters.npoints(2)-1,'observable')`.
- Lines 96-97: computes `fid.sin` using `fid.sin=evolution(spin_system,L_dec,parameters.coil,rho_stack_sin, dw(2),parameters.npoints(2)-1,'observable')`.

### Local helper functions

- Line 102: `grumble()` — `function grumble(spin_system,parameters,H,R,K)`.
  - Representative operation: `if ~ismember(spin_system.bas.formalism,{'sphten-liouv'})`.
  - Representative operation: `error('this function is only available for sphten-liouv formalism.')`.

## Syntax

```matlab
fid=wise(spin_system,parameters,H,R,K)
```

## Parameters / inputs

- parameters.spins -working spins, e.g. {'1H',13C'}
- parameters.hi_pwr -amplitude of high power pulses
- on the high-gamma channel, Hz
- parameters.cp_pwr -amplitude of CP pulse on each
- channel during the CP contact
- time, Hz
- parameters.cp_dur -CP contact time duration, s
- parameters.rho0 -initial state
- parameters.coil -detection state
- parameters.sweep -sweep width, Hz for F1, F2
- parameters.npoints -number of points in F1, F2
- H -Hamiltonian matrix, received from context function
- R -relaxation superoperator, received from context function
- K -kinetics superoperator, received from context function

## Outputs

- fid.sin, fid.cos -sine and cosine components
- of the States quadrature

## Implementation structure

- WISE (WIdeline SEparation) is a powder MAS heteronuclear correlation
- experiment. In the common 1H-13C implementation, molecular dynamics
- information is contained in 1H line shapes that are separated in the
- second dimension by 13C chemical shifts. Further information in:
- fid=wise(spin_system,parameters,H,R,K)
- parameters.spins -working spins, e.g. {'1H',13C'}
- parameters.hi_pwr -amplitude of high power pulses
- on the high-gamma channel, Hz
- parameters.cp_pwr -amplitude of CP pulse on each
- channel during the CP contact
- time, Hz
- parameters.cp_dur -CP contact time duration, s

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `grumble()`, `decouple()`, `operator()`, `speye()`, `step()`, `evolution()`, `ismember()`, `ismatrix()`, `all()`, `isfield()`, `isscalar()`, `isrow()`, `any()`, `iscell()`, `cellfun()`.
