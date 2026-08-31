# kernel/pulses/shaped_pulse_af.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/kernel/pulses/shaped_pulse_af.m`
- Signature: `[rho,traj,P]=shaped_pulse_af(spin_system,L0,Lx,Ly,rho,rf_frq_list,...`
- Total lines: 299

## Purpose

Shaped pulse in amplitude-frequency coordinates using Fokker-Planck formalism (Eqn. 33 in http://dx.doi.org/10.1016/j.jmr.2016.07.005).

## Physical / mathematical content

- Pulse and waveform utilities. These files encode shaped RF pulses, gradient events, rotating-frame transformations, resonator response, and Lie-group integration of time-dependent driven dynamics.
- The file uses a Fokker-Planck-style enlarged state space in which spatial or orientational coordinates are promoted to extra dimensions and coupled to spin dynamics through differential operators.
- Propagation is accelerated with a Krylov-subspace method, replacing direct matrix exponentiation by projection into a much smaller Arnoldi/Lanczos-type subspace.
- Signal processing is central here: the code moves between time and frequency domains, typically using FFT conventions, apodisation, zero filling, or heterodyne frequency shifts.

## Numerical / algorithmic content

- Time propagation is explicit. In Spinach this usually means repeated application of matrix exponentials or propagator factorizations to density operators or state vectors in Hilbert/Liouville/Fokker-Planck space.
- The implementation explicitly addresses performance engineering through parallel or GPU execution, which matters because Spinach operators can become extremely large after basis expansion or powder/spatial lifting.
- A Krylov-subspace or Arnoldi construction is used to avoid forming or exponentiating very large dense propagators directly.
- The file contains an explicit `grumble(...)` validator, which is Spinach convention for front-loading dimension, type, and regime checks before expensive linear-algebra work begins.
- The file also defines local helper function(s): `size()`. This usually means the public entry point is supported by tightly coupled validation or helper logic kept private to the file.

## Code-derived implementation details

### Comment-guided execution stages

- Lines 72-73: Set the defaults; implemented by `if ~exist('method','var'), method='expv'; end`.
- Lines 75-77: Check consistency; implemented by `grumble(L0,Lx,Ly,rho,rf_frq_list,rf_amp_list, rf_dur_list,rf_phi,max_rank,method)`.
- Lines 79-80: Get problem dimensions; implemented by `spc_dim=2*max_rank+1; spn_dim=size(rho,1); stk_dim=size(rho,2)`.
- Lines 86-87: Compute RF phases and Fourier derivative operator; implemented by `[phases,d_dphi]=fourdif(spc_dim,1)`.
- Lines 89-90: Add the overall phase; implemented by `phases=phases+rf_phi`.
- Lines 92-93: Build the background Liouvillian; implemented by `F0=polyadic({{opium(spc_dim,1),L0}})`.
- Lines 95-97: Build RF operators at all phases; implemented by `F1=polyadic({{spdiags(cos(phases),0,spc_dim,spc_dim),Lx}, {spdiags(sin(phases),0,spc_dim,spc_dim),Ly}})`.
- Lines 99-100: Build RF phase turning generator; implemented by `M=polyadic({{d_dphi,opium(spn_dim,1)}})`.
- Lines 102-103: Inflate polyadic representations; implemented by `if ~ismember('polyadic',spin_system.sys.enable)`.
- Lines 109-110: Project the state into the first block; implemented by `proj=zeros(spc_dim,1); proj(1)=1; rho=kron(proj,rho)`.
- Lines 112-113: Upload pertinent things to GPU; implemented by `if ismember('gpu',spin_system.sys.enable)`.
- Lines 121-122: Get the trajectory going; implemented by `if nargout>1`.
- Lines 124-125: Preallocate trajectory array; implemented by `traj=cell(1,numel(rf_dur_list)+1)`.
- Lines 127-128: Fold the first point state stack back into Liouville space; implemented by `traj{1}=squeeze(sum(reshape(full(rho),[spn_dim spc_dim stk_dim]),2))`.
- Lines 132-133: Run the pulse; implemented by `switch method`.
- Lines 137-138: Start feedback timer; implemented by `feedback=tic()`.
- Lines 140-141: Run Krylov propagation; implemented by `for n=1:numel(rf_frq_list)`.
- Lines 143-145: Take a step forward; implemented by `rho=step(spin_system,F0+rf_amp_list(n)*F1+2i*pi*rf_frq_list(n)*M, rho,rf_dur_list(n))`.

### Control flow inferred from the code

- Line 73: conditional branch on `~exist('method','var'), method='expv'; end`.
- Line 103: conditional branch on `~ismember('polyadic',spin_system.sys.enable)`.
- Line 113: conditional branch on `ismember('gpu',spin_system.sys.enable)`.
- Line 122: conditional branch on `nargout>1`.
- Line 133: dispatches on `method`; cases `'expv'`, `'expm'`, `'evolution'`.
- Line 141: `for` loop over `n=1:numel(rf_frq_list)`.
- Line 148: conditional branch on `(n==numel(rf_frq_list))||(toc(feedback)>1)`.
- Line 155: conditional branch on `nargout>1`.
- Line 165: conditional branch on `nargout>2, error('use ''expm'' method to get effective propagator.'); end`.
- Line 170: conditional branch on `nargout>2, P_tot=speye(size(F0)); end`.
- Line 173: `for` loop over `n=1:numel(rf_frq_list)`.
- Line 182: conditional branch on `nargout>1`.
- Line 190: conditional branch on `nargout>2`.
- Line 197: conditional branch on `nargout>2`.

### Key state/data transformations

- Lines 80: computes `spc_dim` using `spc_dim=2*max_rank+1; spn_dim=size(rho,1); stk_dim=size(rho,2)`.
- Lines 87: computes `[phases,d_dphi]` using `[phases,d_dphi]=fourdif(spc_dim,1)`.
- Lines 90: computes `phases` using `phases=phases+rf_phi`.
- Lines 93: computes `F0` using `F0=polyadic({{opium(spc_dim,1),L0}})`.
- Lines 96-97: computes `F1` using `F1=polyadic({{spdiags(cos(phases),0,spc_dim,spc_dim),Lx}, {spdiags(sin(phases),0,spc_dim,spc_dim),Ly}})`.
- Lines 100: computes `M` using `M=polyadic({{d_dphi,opium(spn_dim,1)}})`.
- Lines 110: computes `proj` using `proj=zeros(spc_dim,1); proj(1)=1; rho=kron(proj,rho)`.
- Lines 116: computes `location` using `location='GPU'`.
- Lines 125: computes `traj` using `traj=cell(1,numel(rf_dur_list)+1)`.
- Lines 128: computes `traj{1}` using `traj{1}=squeeze(sum(reshape(full(rho),[spn_dim spc_dim stk_dim]),2))`.
- Lines 138: computes `feedback` using `feedback=tic()`.
- Lines 144-145: computes `rho` using `rho=step(spin_system,F0+rf_amp_list(n)*F1+2i*pi*rf_frq_list(n)*M, rho,rf_dur_list(n))`.
- Lines 158: computes `traj{n+1}` using `traj{n+1}=squeeze(sum(reshape(full(rho),[spn_dim spc_dim stk_dim]),2))`.
- Lines 176: computes `P` using `P=propagator(spin_system,F0+rf_amp_list(n)*F1+2i*pi*rf_frq_list(n)*M,rf_dur_list(n))`.
- Lines 191: computes `P_tot` using `P_tot=clean_up(spin_system,P*P_tot,spin_system.tols.prop_chop)`.
- Lines 249: computes `traj{n}` using `traj{n}=gather(traj{n})`.

### Local helper functions

- Line 258: `grumble()` — `function grumble(L0,Lx,Ly,rho,rf_frq_list,rf_amp_list,`.
  - Representative operation: `rf_dur_list,rf_phi,max_rank,method)`.
  - Representative operation: `if (~isnumeric(L0))||(~isnumeric(Lx))||(~isnumeric(Ly))|| (size(L0,1)~=size(L0,2))||(size(Lx,1)~=size(Lx,2))|| (size(Ly,1)~=size(Ly,2))||(size(L0,1)~=size(Lx,1))|| (size…`.

## Syntax

```matlab
[rho,traj,P]=shaped_pulse_af(spin_system,L0,Lx,Ly,rho,rf_frq_list,...
rf_amp_list,rf_dur_list,rf_phi,...
max_rank,method)
```

## Parameters / inputs

- L0 -drift Liouvillian that continues
- running in the background
- Lx -X projection of the RF operator
- Ly -Y projection of the RF operator
- rho -initial state vector or a horizontal
- stack thereof
- rf_frq_list -a vector of RF frequencies at each
- time slice (relative to the offsets
- and/or rotating frames that were
- used to make the background L0 that
- you have supplied, Hz
- rf_amp_list -a vector of RF amplitudes at each
- time slice, rad/s
- rf_dur_list -a vector of time slice durations,
- in seconds
- rf_phi -RF phase of the first pulse slice
- max_rank -maximum rank of the Fokker-Planck
- theory, increase until the answer
- stops changing, 2 is a good start
- method -propagation method, 'expv' for Krylov
- propagation, 'expm' for exponential
- propagation, 'evolution' for Spinach
- evolution function

## Outputs

- rho -final state vector or a stack thereof
- traj -system trajectory as a [1 x (nsteps+1)]
- cell array, the first point is the ini-
- tial condition
- P -effective pulse propagator (expensive),
- only available for the 'expm' method
- Note: the pulse is assumed to be piecewise-constant and should be
- supplied with sufficiently fine time discretisation to pro-
- perly reproduce the waveform.
- Note: make it dead certain that your freqiency has the correct
- sign; wrong sign means that the pulse hits very far away
- from your intended location. This is the principal source
- of bugs when using this function.

## Implementation structure

- Shaped pulse in amplitude-frequency coordinates using Fokker-Planck
- formalism (Eqn. 33 in http://dx.doi.org/10.1016/j.jmr.2016.07.005).
- [rho,traj,P]=shaped_pulse_af(spin_system,L0,Lx,Ly,rho,rf_frq_list,...
- rf_amp_list,rf_dur_list,rf_phi,...
- max_rank,method)
- L0 -drift Liouvillian that continues
- running in the background
- Lx -X projection of the RF operator
- Ly -Y projection of the RF operator
- rho -initial state vector or a horizontal
- stack thereof
- rf_frq_list -a vector of RF frequencies at each

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `exist()`, `grumble()`, `report()`, `num2str()`, `fourdif()`, `polyadic()`, `opium()`, `spdiags()`, `ismember()`, `complex()`, `inflate()`, `proj()`, `gpuArray()`, `squeeze()`, `tic()`, `step()`.
