# experiments/imaging/udd_dec.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/experiments/imaging/udd_dec.m`
- Signature: `mri=udd_dec(spin_system,parameters,H,R,K,G,F)`
- Total lines: 118

## Purpose

The effect of Uhrid Dynamic Decoupling (UDD) pulse sequence on the MRI phantom. The function runs the UDD and then pro- jects out the user-specified spin state, returning the cor- responding image. Syntax: mri=udd_dec(spin_system,parameters,H,R,K,G,F) This sequence must be called from the imaging() context, which would provide H,R,K,G, and F. Parameters: parameters.dec_time -total duration of the sequence parameters.

## Physical / mathematical content

- Imaging sequence implementations. They build spatially resolved Liouvillians that include gradients, slice-selection RF terms, diffusion, and acquisition operators.

## Numerical / algorithmic content

- Time propagation is explicit. In Spinach this usually means repeated application of matrix exponentials or propagator factorizations to density operators or state vectors in Hilbert/Liouville/Fokker-Planck space.
- The file contains an explicit `grumble(...)` validator, which is Spinach convention for front-loading dimension, type, and regime checks before expensive linear-algebra work begins.
- The file also defines local helper function(s): `grumble()`. This usually means the public entry point is supported by tightly coupled validation or helper logic kept private to the file.

## Code-derived implementation details

### Comment-guided execution stages

- Lines 34-35: Check consistency; implemented by `grumble(spin_system,parameters,H,R,K,G,F)`.
- Lines 37-38: Assemble the background; implemented by `B=H+F+1i*R+1i*K`.
- Lines 40-41: Make pulse operators; implemented by `Lp=operator(spin_system,'L+',parameters.spins{1})`.
- Lines 45-46: Apply 90-degree pulse; implemented by `rho=step(spin_system,Ly,parameters.rho0,pi/2)`.
- Lines 48-49: Get UDD delays; implemented by `udd_delays=uhrig_times(parameters.dec_time,parameters.npulses)`.
- Lines 51-52: Run UDD echo train; implemented by `for n=1:(numel(udd_delays)-1)`.
- Lines 54-55: Apply the delay; implemented by `rho=step(spin_system,B,rho,udd_delays(n))`.
- Lines 57-58: Apply the pulse; implemented by `rho=step(spin_system,Lx,rho,pi)`.
- Lines 60-62: Inform the user; implemented by `report(spin_system,['UDD pi pulse ' num2str(n) ' out of ' num2str(numel(udd_delays)-1) ' done.'])`.
- Lines 66-67: Apply the last delay; implemented by `rho=step(spin_system,B,rho,udd_delays(end))`.
- Lines 69-70: Project out the spin state and get its phantom; implemented by `mri=fpl2phan(rho,parameters.coil_st{1},parameters.npts)`.

### Control flow inferred from the code

- Line 52: `for` loop over `n=1:(numel(udd_delays)-1)`.

### Key state/data transformations

- Lines 38: computes `B` using `B=H+F+1i*R+1i*K`.
- Lines 41: computes `Lp` using `Lp=operator(spin_system,'L+',parameters.spins{1})`.
- Lines 42: computes `Ly` using `Ly=kron(speye(prod(parameters.npts)),(Lp-Lp')/2i)`.
- Lines 43: computes `Lx` using `Lx=kron(speye(prod(parameters.npts)),(Lp+Lp')/2)`.
- Lines 46: computes `rho` using `rho=step(spin_system,Ly,parameters.rho0,pi/2)`.
- Lines 49: computes `udd_delays` using `udd_delays=uhrig_times(parameters.dec_time,parameters.npulses)`.
- Lines 70: computes `mri` using `mri=fpl2phan(rho,parameters.coil_st{1},parameters.npts)`.

### Local helper functions

- Line 75: `grumble()` — `function grumble(spin_system,parameters,H,R,K,G,F)`.
  - Representative operation: `if ~ismember(spin_system.bas.formalism,{'sphten-liouv','zeeman-liouv'})`.
  - Representative operation: `error('this function is only available in sphten-liouv and zeeman-liouv formalisms.')`.

## Outputs

- mri -amplitude of the detection state at each point of the
- sample
- Note: the spin state to be observed should be specified in
- parameters.coil_st, the coil phantom is ignored.

## Implementation structure

- The effect of Uhrid Dynamic Decoupling (UDD) pulse sequence
- on the MRI phantom. The function runs the UDD and then pro-
- jects out the user-specified spin state, returning the cor-
- responding image. Syntax:
- mri=udd_dec(spin_system,parameters,H,R,K,G,F)
- This sequence must be called from the imaging() context, which
- would provide H,R,K,G, and F. Parameters:
- parameters.dec_time -total duration of the sequence
- parameters.npulses -number of pulses in the sequence,
- excluding the first pi/2 pulse
- parameters.spins -nuclei on which the sequence
- is to act, e.g. {'1H'}

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `grumble()`, `operator()`, `speye()`, `step()`, `uhrig_times()`, `udd_delays()`, `report()`, `num2str()`, `fpl2phan()`, `ismember()`, `ismatrix()`, `all()`, `iscell()`, `isfield()`, `isscalar()`.
