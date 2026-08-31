# experiments/imaging/cpmg_dec.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/experiments/imaging/cpmg_dec.m`
- Signature: `mri=cpmg_dec(spin_system,parameters,H,R,K,G,F)`
- Total lines: 145

## Purpose

The effect of Carr-Purcell-Meiboom-Gill (CPMG) pulse sequence on the MRI phantom. The function runs the CPMG and then pro- jects out the user-specified spin state, returning the corres- ponding image. Syntax: mri=cpmg_dec(spin_system,parameters,H,R,K,G,F) This sequence must be called from the imaging() context, which would provide H,R,K,G, and F. Parameters: parameters.dec_time -total duration of the sequence paramet

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
- Lines 45-46: Get CPMG delay; implemented by `cpmg_delay=parameters.dec_time/parameters.npulses`.
- Lines 48-49: Apply the first pulse; implemented by `rho=step(spin_system,Ly,parameters.rho0,pi/2)`.
- Lines 51-52: Apply the first delay; implemented by `rho=step(spin_system,B,rho,cpmg_delay/2)`.
- Lines 54-55: Run CPMG echo train; implemented by `for n=1:(parameters.npulses-1)`.
- Lines 57-58: Apply the pulse; implemented by `rho=step(spin_system,Lx,rho,pi)`.
- Lines 60-61: Apply the delay; implemented by `rho=step(spin_system,B,rho,cpmg_delay)`.
- Lines 63-65: Inform the user; implemented by `report(spin_system,['CPMG pi pulse ' num2str(n) ' out of ' num2str(parameters.npulses) ' done.'])`.
- Lines 69-70: Apply the last pulse; implemented by `rho=step(spin_system,Lx,rho,pi)`.
- Lines 72-73: Apply the last delay; implemented by `rho=step(spin_system,B,rho,cpmg_delay/2)`.
- Lines 75-76: Project out the spin state and get its phantom; implemented by `mri=fpl2phan(rho,parameters.coil_st{1},parameters.npts)`.

### Control flow inferred from the code

- Line 55: `for` loop over `n=1:(parameters.npulses-1)`.

### Key state/data transformations

- Lines 38: computes `B` using `B=H+F+1i*R+1i*K`.
- Lines 41: computes `Lp` using `Lp=operator(spin_system,'L+',parameters.spins{1})`.
- Lines 42: computes `Ly` using `Ly=kron(speye(prod(parameters.npts)),(Lp-Lp')/2i)`.
- Lines 43: computes `Lx` using `Lx=kron(speye(prod(parameters.npts)),(Lp+Lp')/2)`.
- Lines 46: computes `cpmg_delay` using `cpmg_delay=parameters.dec_time/parameters.npulses`.
- Lines 49: computes `rho` using `rho=step(spin_system,Ly,parameters.rho0,pi/2)`.
- Lines 76: computes `mri` using `mri=fpl2phan(rho,parameters.coil_st{1},parameters.npts)`.

### Local helper functions

- Line 81: `grumble()` — `function grumble(spin_system,parameters,H,R,K,G,F)`.
  - Representative operation: `if ~ismember(spin_system.bas.formalism,{'sphten-liouv','zeeman-liouv'})`.
  - Representative operation: `error('this function is only available in sphten-liouv and zeeman-liouv formalisms.')`.

## Outputs

- mri -amplitude of the detection state at each point of the
- sample
- Note: the spin state to be observed should be specified in
- parameters.coil_st, the coil phantom is ignored.

## Implementation structure

- The effect of Carr-Purcell-Meiboom-Gill (CPMG) pulse sequence
- on the MRI phantom. The function runs the CPMG and then pro-
- jects out the user-specified spin state, returning the corres-
- ponding image. Syntax:
- mri=cpmg_dec(spin_system,parameters,H,R,K,G,F)
- This sequence must be called from the imaging() context, which
- would provide H,R,K,G, and F. Parameters:
- parameters.dec_time -total duration of the sequence
- parameters.npulses -number of pulses in the sequence,
- excluding the first pi/2 pulse
- parameters.spins -nuclei on which the sequence
- is to act, e.g. {'1H'}

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `grumble()`, `operator()`, `speye()`, `step()`, `report()`, `num2str()`, `fpl2phan()`, `ismember()`, `ismatrix()`, `all()`, `iscell()`, `isfield()`, `ischar()`, `isvector()`, `any()`, `isscalar()`.
