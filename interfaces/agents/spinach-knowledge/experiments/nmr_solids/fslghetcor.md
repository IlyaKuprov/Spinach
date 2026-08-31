# experiments/nmr_solids/fslghetcor.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/experiments/nmr_solids/fslghetcor.m`
- Signature: `fid=fslghetcor(spin_system,parameters,H,R,K)`
- Total lines: 235

## Purpose

Heteronuclear correlation MAS NMR experiment with frequency-switched Lee-Goldburg homonuclear decoupling. Further details in:

## Physical / mathematical content

- Solid-state pulse sequence implementations. The core ingredients are anisotropic Hamiltonians, rotor synchronisation, cross-polarisation, recoupling/decoupling, and powder or rotor-stack propagation.
- The file uses a Fokker-Planck-style enlarged state space in which spatial or orientational coordinates are promoted to extra dimensions and coupled to spin dynamics through differential operators.
- Orientation or trajectory averaging is performed numerically, so grid design, weights, and integration error control matter directly to accuracy and runtime.

## Numerical / algorithmic content

- Time propagation is explicit. In Spinach this usually means repeated application of matrix exponentials or propagator factorizations to density operators or state vectors in Hilbert/Liouville/Fokker-Planck space.
- The implementation explicitly addresses performance engineering through parallel or GPU execution, which matters because Spinach operators can become extremely large after basis expansion or powder/spatial lifting.
- Numerical integration over angles or geometry is part of the implementation, so point placement and weights are as important as the local Hamiltonian calculations.
- The file contains an explicit `grumble(...)` validator, which is Spinach convention for front-loading dimension, type, and regime checks before expensive linear-algebra work begins.
- The file also defines local helper function(s): `grumble()`. This usually means the public entry point is supported by tightly coupled validation or helper logic kept private to the file.

## Code-derived implementation details

### Comment-guided execution stages

- Lines 58-59: Consistency enforcement; implemented by `grumble(spin_system,parameters,H,R,K)`.
- Lines 61-62: Build 1H and 13C control operators; implemented by `Hx=operator(spin_system,'Lx',parameters.spins{1}); Hx=kron(speye(parameters.spc_dim),Hx)`.
- Lines 67-68: Liouvillian; implemented by `L=H+1i*R+1i*K`.
- Lines 70-71: Proton pulse evolution generators; implemented by `L_HPx=L+2*pi*parameters.hi_pwr*Hx; L_HPy=L+2*pi*parameters.hi_pwr*Hy`.
- Lines 74-76: FSLG pulse evolution generators; implemented by `L_FSLG_p=L-2*pi*parameters.offset(1)*Hz +2*pi*(parameters.hi_pwr/sqrt(2))*Hz`.
- Lines 84-86: CP period evolution generator; implemented by `L_c=L-2*pi*parameters.cp_pwr(1)*Hy +2*pi*parameters.cp_pwr(2)*Cx`.
- Lines 88-89: Sequence timing shorthands; implemented by `hi_pwr_90deg=1/(4*parameters.hi_pwr)`.
- Lines 93-94: Flip by 90 degrees + magic angle; implemented by `rho_cos=step(spin_system,L_HPx,parameters.rho0,hi_pwr_90deg+hi_pwr_magic)`.
- Lines 97-98: Preallocate and start the F1 trajectories; implemented by `traj_cos=zeros(numel(rho_cos),parameters.npoints(1),'like',1i)`.
- Lines 102-103: Compute the F1 trajectory, cos part; implemented by `if ismember('gpu',spin_system.sys.enable)`.
- Lines 115-116: Compute the F1 trajectory, sin part; implemented by `if ismember('gpu',spin_system.sys.enable)`.
- Lines 128-129: Flip the 1H back from the magic angle to the x,y-plane; implemented by `traj_cos=step(spin_system,L_HMx,traj_cos,hi_pwr_magic)`.
- Lines 132-133: Save work by stacking cos and sin parts; implemented by `traj_cos_sin=[traj_cos traj_sin]`.
- Lines 135-136: Run the CP contact period; implemented by `traj_cos_sin=step(spin_system,L_c,traj_cos_sin,parameters.cp_dur)`.
- Lines 138-139: Decouple 1H for the acquisition period; implemented by `[L_dec,traj_cos_sin]=decouple(spin_system,L,traj_cos_sin,parameters.spins(1))`.
- Lines 141-143: Run the F2 evolution; implemented by `fid_cos_sin=evolution(spin_system,L_dec,parameters.coil,traj_cos_sin, dwell_times(2),parameters.npoints(2)-1,'observable')`.
- Lines 145-146: Unstack cos and sin parts; implemented by `fid.cos=fid_cos_sin(:,1:(end/2))`.

### Control flow inferred from the code

- Line 103: conditional branch on `ismember('gpu',spin_system.sys.enable)`.
- Line 106: `for` loop over `k=2:parameters.npoints(1)`.
- Line 107: `for` loop over `n=1:parameters.nblocks`.
- Line 116: conditional branch on `ismember('gpu',spin_system.sys.enable)`.
- Line 119: `for` loop over `k=2:parameters.npoints(1)`.
- Line 120: `for` loop over `n=1:parameters.nblocks`.

### Key state/data transformations

- Lines 62: computes `Hx` using `Hx=operator(spin_system,'Lx',parameters.spins{1}); Hx=kron(speye(parameters.spc_dim),Hx)`.
- Lines 63: computes `Hy` using `Hy=operator(spin_system,'Ly',parameters.spins{1}); Hy=kron(speye(parameters.spc_dim),Hy)`.
- Lines 64: computes `Hz` using `Hz=operator(spin_system,'Lz',parameters.spins{1}); Hz=kron(speye(parameters.spc_dim),Hz)`.
- Lines 65: computes `Cx` using `Cx=operator(spin_system,'Lx',parameters.spins{2}); Cx=kron(speye(parameters.spc_dim),Cx)`.
- Lines 68: computes `L` using `L=H+1i*R+1i*K`.
- Lines 71: computes `L_HPx` using `L_HPx=L+2*pi*parameters.hi_pwr*Hx; L_HPy=L+2*pi*parameters.hi_pwr*Hy`.
- Lines 72: computes `L_HMx` using `L_HMx=L-2*pi*parameters.hi_pwr*Hx; L_HMy=L-2*pi*parameters.hi_pwr*Hy`.
- Lines 75-76: computes `L_FSLG_p` using `L_FSLG_p=L-2*pi*parameters.offset(1)*Hz +2*pi*(parameters.hi_pwr/sqrt(2))*Hz`.
- Lines 77-78: computes `L_FSLG_m` using `L_FSLG_m=L-2*pi*parameters.offset(1)*Hz -2*pi*(parameters.hi_pwr/sqrt(2))*Hz`.
- Lines 79: computes `L1_cos` using `L1_cos=L_FSLG_m+2*pi*parameters.hi_pwr*Hy`.
- Lines 80: computes `L2_cos` using `L2_cos=L_FSLG_p-2*pi*parameters.hi_pwr*Hy`.
- Lines 81: computes `L1_sin` using `L1_sin=L_FSLG_p+2*pi*parameters.hi_pwr*Hx`.
- Lines 82: computes `L2_sin` using `L2_sin=L_FSLG_m-2*pi*parameters.hi_pwr*Hx`.
- Lines 85-86: computes `L_c` using `L_c=L-2*pi*parameters.cp_pwr(1)*Hy +2*pi*parameters.cp_pwr(2)*Cx`.
- Lines 89: computes `hi_pwr_90deg` using `hi_pwr_90deg=1/(4*parameters.hi_pwr)`.
- Lines 90: computes `hi_pwr_magic` using `hi_pwr_magic=acos(1/sqrt(3))/(2*pi)/parameters.hi_pwr`.
- Lines 91: computes `dwell_times` using `dwell_times=1./parameters.sweep`.
- Lines 94: computes `rho_cos` using `rho_cos=step(spin_system,L_HPx,parameters.rho0,hi_pwr_90deg+hi_pwr_magic)`.

### Local helper functions

- Line 152: `grumble()` — `function grumble(spin_system,parameters,H,R,K)`.
  - Representative operation: `if ~ismember(spin_system.bas.formalism,{'sphten-liouv'})`.
  - Representative operation: `error('this function is only available for sphten-liouv formalism.')`.

## Syntax

```matlab
fid=fslghetcor(spin_system,parameters,H,R,K)
```

## Parameters / inputs

- parameters.spins -working spins, e.g. {'1H',13C'}
- parameters.hi_pwr -amplitude of high power pulses
- on the high-gamma channel, Hz
- parameters.cp_pwr -amplitude of CP pulse on each
- channel during the CP contact
- time, Hz
- parameters.cp_dur -CP contact time duration, s
- parameters.offset -transmitter offsets on the
- two channels, Hz
- parameters.nblocks -number of FSLG blocks per
- indirect-dimension point
- parameters.spc_dim -Fokker-Planck spatial dimension
- parameters.rho0 -initial state
- parameters.coil -detection state
- parameters.sweep -sweep width, Hz for F1, F2
- parameters.npoints -number of points in F1, F2
- H -Hamiltonian superoperator, received from context function
- R -relaxation superoperator, received from context function
- K -kinetics superoperator, received from context function

## Outputs

- fid.sin, fid.cos -sine and cosine components
- of the States quadrature

## Implementation structure

- Heteronuclear correlation MAS NMR experiment with frequency-switched
- Lee-Goldburg homonuclear decoupling. Further details in:
- fid=fslghetcor(spin_system,parameters,H,R,K)
- parameters.spins -working spins, e.g. {'1H',13C'}
- parameters.hi_pwr -amplitude of high power pulses
- on the high-gamma channel, Hz
- parameters.cp_pwr -amplitude of CP pulse on each
- channel during the CP contact
- time, Hz
- parameters.cp_dur -CP contact time duration, s
- parameters.offset -transmitter offsets on the
- two channels, Hz

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `grumble()`, `operator()`, `speye()`, `acos()`, `step()`, `traj_cos()`, `traj_sin()`, `ismember()`, `gpuArray()`, `clear()`, `decouple()`, `evolution()`, `dwell_times()`, `fid_cos_sin()`, `ismatrix()`, `all()`.
