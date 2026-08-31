# experiments/nmr_solids/cn2d_dq.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/experiments/nmr_solids/cn2d_dq.m`
- Signature: `fid=cn2d_dq(spin_system,parameters,H,R,K)`
- Total lines: 141

## Purpose

Double-quantum version of the 13C-detected 14N-13C MAS 2D correlation experiment described by Jarvis, Haies, Williamson and Carravetta in

## Physical / mathematical content

- Solid-state pulse sequence implementations. The core ingredients are anisotropic Hamiltonians, rotor synchronisation, cross-polarisation, recoupling/decoupling, and powder or rotor-stack propagation.
- The file uses a Fokker-Planck-style enlarged state space in which spatial or orientational coordinates are promoted to extra dimensions and coupled to spin dynamics through differential operators.
- Orientation or trajectory averaging is performed numerically, so grid design, weights, and integration error control matter directly to accuracy and runtime.

## Numerical / algorithmic content

- Time propagation is explicit. In Spinach this usually means repeated application of matrix exponentials or propagator factorizations to density operators or state vectors in Hilbert/Liouville/Fokker-Planck space.
- Numerical integration over angles or geometry is part of the implementation, so point placement and weights are as important as the local Hamiltonian calculations.
- The file contains an explicit `grumble(...)` validator, which is Spinach convention for front-loading dimension, type, and regime checks before expensive linear-algebra work begins.
- The file also defines local helper function(s): `grumble()`. This usually means the public entry point is supported by tightly coupled validation or helper logic kept private to the file.

## Code-derived implementation details

### Comment-guided execution stages

- Lines 45-46: Check consistency; implemented by `grumble(spin_system,parameters,H,R,K)`.
- Lines 48-49: Compose the Liouvillian; implemented by `L=H+1i*R+1i*K`.
- Lines 51-52: Get evolution timesteps; implemented by `timestep=1./parameters.sweep`.
- Lines 54-55: Get pulse operators; implemented by `Cp=operator(spin_system,'L+',parameters.spins{2})`.
- Lines 61-62: Apply 14N pulse; implemented by `rho_cos=evolution(spin_system,L+2*pi*parameters.rf_pwr*Nx,[],parameters.rho0,parameters.rf_dur,1,'final')`.
- Lines 65-66: Apply coherence selection; implemented by `rho_cos=coherence(spin_system,rho_cos,{{parameters.spins{2},[-1 1]},{parameters.spins{1},[-2 2]}})`.
- Lines 69-70: Run F1 evolution and 13C pulse; implemented by `rho_stack_cos=evolution(spin_system,L,[],rho_cos,timestep(1)/2,parameters.npoints(1)-1,'trajectory')`.
- Lines 76-77: Apply 14N pulse; implemented by `rho_stack_cos=evolution(spin_system,L+2*pi*parameters.rf_pwr*Nx,[],rho_stack_cos,parameters.rf_dur,1,'final')`.
- Lines 80-81: Run F2 evolution; implemented by `fid.cos=evolution(spin_system,L,parameters.coil,rho_stack_cos,timestep(2),parameters.npoints(2)-1,'observable')`.

### Key state/data transformations

- Lines 49: computes `L` using `L=H+1i*R+1i*K`.
- Lines 52: computes `timestep` using `timestep=1./parameters.sweep`.
- Lines 55: computes `Cp` using `Cp=operator(spin_system,'L+',parameters.spins{2})`.
- Lines 56: computes `Np` using `Np=operator(spin_system,'L+',parameters.spins{1})`.
- Lines 57: computes `Cx` using `Cx=(Cp+Cp')/2; Cx=kron(speye(parameters.spc_dim),Cx)`.
- Lines 58: computes `Nx` using `Nx=(Np+Np')/2; Nx=kron(speye(parameters.spc_dim),Nx)`.
- Lines 59: computes `Ny` using `Ny=(Np-Np')/2i; Ny=kron(speye(parameters.spc_dim),Ny)`.
- Lines 62: computes `rho_cos` using `rho_cos=evolution(spin_system,L+2*pi*parameters.rf_pwr*Nx,[],parameters.rho0,parameters.rf_dur,1,'final')`.
- Lines 63: computes `rho_sin` using `rho_sin=evolution(spin_system,L+2*pi*parameters.rf_pwr*(Nx+Ny)/sqrt(2),[],parameters.rho0,parameters.rf_dur,1,'final')`.
- Lines 70: computes `rho_stack_cos` using `rho_stack_cos=evolution(spin_system,L,[],rho_cos,timestep(1)/2,parameters.npoints(1)-1,'trajectory')`.
- Lines 71: computes `rho_stack_sin` using `rho_stack_sin=evolution(spin_system,L,[],rho_sin,timestep(1)/2,parameters.npoints(1)-1,'trajectory')`.
- Lines 81: computes `fid.cos` using `fid.cos=evolution(spin_system,L,parameters.coil,rho_stack_cos,timestep(2),parameters.npoints(2)-1,'observable')`.
- Lines 82: computes `fid.sin` using `fid.sin=evolution(spin_system,L,parameters.coil,rho_stack_sin,timestep(2),parameters.npoints(2)-1,'observable')`.

### Local helper functions

- Line 87: `grumble()` — `function grumble(spin_system,parameters,H,R,K)`.
  - Representative operation: `if ~ismember(spin_system.bas.formalism,{'sphten-liouv'})`.
  - Representative operation: `error('this function is only available for sphten-liouv formalism.')`.

## Syntax

```matlab
fid=cn2d_dq(spin_system,parameters,H,R,K)
```

## Parameters / inputs

- parameters.spins -isotopes to which the sequence is
- applied, specified as a cell array
- with 14N first, and 13C second
- parameters.spc_dim -Fokker-Planck spatial dimension
- parameters.sweep -sweep widths in the two dimensions, Hz
- parameters.npoints -numbers of points in the two dimensions
- parameters.rho0 -initial state
- parameters.coil -detection state
- parameters.rf_pwr -RF power on 14N, Hz
- parameters.rf_dur -RF pulse duration on 14N, seconds

## Outputs

- fid.sin
- fid.cos -sine and cosine components
- of the States quadrature

## Implementation structure

- Double-quantum version of the 13C-detected 14N-13C MAS 2D correlation
- experiment described by Jarvis, Haies, Williamson and Carravetta in
- fid=cn2d_dq(spin_system,parameters,H,R,K)
- parameters.spins -isotopes to which the sequence is
- applied, specified as a cell array
- with 14N first, and 13C second
- parameters.spc_dim -Fokker-Planck spatial dimension
- parameters.sweep -sweep widths in the two dimensions, Hz
- parameters.npoints -numbers of points in the two dimensions
- parameters.rho0 -initial state
- parameters.coil -detection state
- parameters.rf_pwr -RF power on 14N, Hz

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `grumble()`, `operator()`, `speye()`, `evolution()`, `coherence()`, `timestep()`, `step()`, `ismember()`, `ismatrix()`, `all()`, `isfield()`, `elseif()`.
