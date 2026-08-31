# experiments/spen/st_ideal.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/experiments/spen/st_ideal.m`
- Signature: `inten=st_ideal(spin_system,parameters,H,R,K,G,F)`
- Total lines: 132

## Purpose

The ideal Stejskal-Tanner pulse sequence using the notation from Figure 1 in http://dx.doi.org/0.1002/cmr.a.21241 with no gaps be- tween pulse sequence events. Syntax: inten=st_ideal(spin_system,parameters,H,R,K,G,F) This sequence must be called from the imaging() context, which would provide H, R, K, G, and F. Parameters: parameters.spins -working spin. parameters.g_amp -gradient amplitude, T/m parameters.delta_sml 

## Physical / mathematical content

- SPEN experiment implementations. These files combine shaped pulses, gradients, spatial encoding, and often diffusion-aware propagation.

## Numerical / algorithmic content

- Time propagation is explicit. In Spinach this usually means repeated application of matrix exponentials or propagator factorizations to density operators or state vectors in Hilbert/Liouville/Fokker-Planck space.
- The file contains an explicit `grumble(...)` validator, which is Spinach convention for front-loading dimension, type, and regime checks before expensive linear-algebra work begins.
- The file also defines local helper function(s): `grumble()`. This usually means the public entry point is supported by tightly coupled validation or helper logic kept private to the file.

## Code-derived implementation details

### Comment-guided execution stages

- Lines 35-36: Check consistency; implemented by `grumble(spin_system,parameters,H,R,K,F,G)`.
- Lines 38-39: Compose Liouvillian; implemented by `L=H+F+1i*R+1i*K`.
- Lines 41-42: Get pulse operators; implemented by `Lp=operator(spin_system,'L+',parameters.spins{1})`.
- Lines 45-46: Apply the excitation pulse; implemented by `rho=step(spin_system,Ly,parameters.rho0,pi/2)`.
- Lines 48-49: Evolve under the first gradient; implemented by `rho=evolution(spin_system,L+parameters.g_amp*G{1},[],rho,parameters.delta_sml,1,'final')`.
- Lines 51-52: Diffusion delay; implemented by `rho=step(spin_system,L,rho,(parameters.delta_big-parameters.delta_sml)/2)`.
- Lines 54-55: Apply the refocusing pulse; implemented by `rho=step(spin_system,Ly,rho,pi)`.
- Lines 57-58: Evolve to the second gradient; implemented by `rho=step(spin_system,L,rho,(parameters.delta_big-parameters.delta_sml)/2)`.
- Lines 60-61: Evolve under the second gradient; implemented by `rho=evolution(spin_system,L+parameters.g_amp*G{1},[],rho,parameters.delta_sml,1,'final')`.
- Lines 63-64: Intensity is the first point in the FID; implemented by `inten=abs(parameters.coil'*rho)`.

### Key state/data transformations

- Lines 39: computes `L` using `L=H+F+1i*R+1i*K`.
- Lines 42: computes `Lp` using `Lp=operator(spin_system,'L+',parameters.spins{1})`.
- Lines 43: computes `Ly` using `Ly=kron(speye(prod(parameters.npts)),(Lp-Lp')/2i)`.
- Lines 46: computes `rho` using `rho=step(spin_system,Ly,parameters.rho0,pi/2)`.
- Lines 64: computes `inten` using `inten=abs(parameters.coil'*rho)`.

### Local helper functions

- Line 69: `grumble()` — `function grumble(spin_system,parameters,H,R,K,F,G)`.
  - Representative operation: `if ~ismember(spin_system.bas.formalism,{'sphten-liouv'})`.
  - Representative operation: `error('this function is only available for sphten-liouv formalism.')`.

## Outputs

- inten -the absolute value of the first point in
- the free induction decay; this number is
- proportional to the integral of the real
- part of the correctly phased spectrum

## Implementation structure

- The ideal Stejskal-Tanner pulse sequence using the notation from
- Figure 1 in http://dx.doi.org/0.1002/cmr.a.21241 with no gaps be-
- tween pulse sequence events. Syntax:
- inten=st_ideal(spin_system,parameters,H,R,K,G,F)
- This sequence must be called from the imaging() context, which
- would provide H, R, K, G, and F. Parameters:
- parameters.spins -working spin.
- parameters.g_amp -gradient amplitude, T/m
- parameters.delta_sml -the small delta parameter
- (see the figure)
- parameters.delta_big -the big delta parameter
- inten -the absolute value of the first point in

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `grumble()`, `operator()`, `speye()`, `step()`, `evolution()`, `ismember()`, `ismatrix()`, `all()`, `iscell()`, `isfield()`, `elseif()`, `isfloat()`, `eps()`.
