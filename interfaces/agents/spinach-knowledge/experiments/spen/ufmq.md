# experiments/spen/ufmq.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/experiments/spen/ufmq.m`
- Signature: `fid=ufmq(spin_system,parameters,H,R,K,G,F)`
- Total lines: 262

## Purpose

Ultrafast multiple-quantum NMR, a literal implementation of Figure 1A from (http://dx.doi.org/10.1002/cphc.201800667). Syntax: fid=ufmq_nmr(spin_system,parameters,H,R,K,G,F) This sequence must be called from the imaging() context, which would provide H, R, K, G, and F. Parameters: parameters.spins nuclei on which the sequence runs parameters.dims size of the sample, m parameters.npts number of grid points parameters.

## Physical / mathematical content

- SPEN experiment implementations. These files combine shaped pulses, gradients, spatial encoding, and often diffusion-aware propagation.

## Numerical / algorithmic content

- Time propagation is explicit. In Spinach this usually means repeated application of matrix exponentials or propagator factorizations to density operators or state vectors in Hilbert/Liouville/Fokker-Planck space.
- The implementation explicitly addresses performance engineering through parallel or GPU execution, which matters because Spinach operators can become extremely large after basis expansion or powder/spatial lifting.
- The file contains an explicit `grumble(...)` validator, which is Spinach convention for front-loading dimension, type, and regime checks before expensive linear-algebra work begins.
- The file also defines local helper function(s): `grumble()`. This usually means the public entry point is supported by tightly coupled validation or helper logic kept private to the file.

## Code-derived implementation details

### Comment-guided execution stages

- Lines 51-52: Check consistency; implemented by `grumble(spin_system,parameters,H,R,K,G,F)`.
- Lines 54-55: Compose Liouvillian; implemented by `L=H+F+1i*R+1i*K`.
- Lines 57-58: Get pulse operators; implemented by `Lp=operator(spin_system,'L+',parameters.spins{1})`.
- Lines 65-67: Get chirp pulse waveform; implemented by `[Cx,Cy]=chirp_pulse(parameters.pulsenpoints,parameters.Te, parameters.BW,parameters.nWURST,parameters.chirptype)`.
- Lines 70-71: Get the gradient pulse waveform for spatial encoding; implemented by `gradient_amplitudes=parameters.Ge(1)*ones(parameters.pulsenpoints,1)`.
- Lines 73-74: Apply the first pulse; implemented by `rho=step(spin_system,Lx,parameters.rho0,pi/2)`.
- Lines 76-77: Apply first delay; implemented by `rho=evolution(spin_system,L,[],rho,parameters.delay,1,'final')`.
- Lines 79-80: Apply 180 deg pulse; implemented by `rho=step(spin_system,Lx,rho,pi)`.
- Lines 82-83: Apply second delay; implemented by `rho=evolution(spin_system,L,[],rho,parameters.delay,1,'final')`.
- Lines 85-86: Select operator for the second 90 deg pulse; implemented by `if mod(parameters.mqorder,2)==0`.
- Lines 88-89: Apply the second 90 deg pulse on x; implemented by `rho=step(spin_system,Lx,rho,pi/2)`.
- Lines 93-94: Apply the second 90 deg pulse on y; implemented by `rho=step(spin_system,Ly,rho,pi/2)`.
- Lines 98-99: Apply the pair of chirp pulses with opposite gradients; implemented by `rho=coherence(spin_system,rho,{{parameters.spins{1},+parameters.mqorder}})`.
- Lines 107-108: Apply third 90 deg pulse; implemented by `rho=step(spin_system,Lx,rho,pi/2)`.
- Lines 110-111: Coherence selection; implemented by `rho=coherence(spin_system,rho,{{parameters.spins{1},+1}})`.
- Lines 113-114: Apply prephasing (to set the echo in the right position); implemented by `Taq=parameters.npoints*parameters.deltat`.
- Lines 117-118: Preallocate the fid; implemented by `fid=zeros([parameters.npoints parameters.nloops],'like',1i)`.
- Lines 120-121: Upload to GPU; implemented by `if ismember('gpu',spin_system.sys.enable)`.

### Control flow inferred from the code

- Line 61: conditional branch on `~ismember('polyadic',spin_system.sys.enable)`.
- Line 86: conditional branch on `mod(parameters.mqorder,2)==0`.
- Line 121: conditional branch on `ismember('gpu',spin_system.sys.enable)`.
- Line 131: `for` loop over `m=1:parameters.nloops`.
- Line 140: `for` loop over `k=1:parameters.npoints`.
- Line 153: conditional branch on `ismember('gpu',spin_system.sys.enable)`.

### Key state/data transformations

- Lines 55: computes `L` using `L=H+F+1i*R+1i*K`.
- Lines 58: computes `Lp` using `Lp=operator(spin_system,'L+',parameters.spins{1})`.
- Lines 59: computes `Lx` using `Lx=polyadic({{speye(prod(parameters.npts)),(Lp+Lp')/2}})`.
- Lines 60: computes `Ly` using `Ly=polyadic({{speye(prod(parameters.npts)),(Lp-Lp')/2i}})`.
- Lines 66-67: computes `[Cx,Cy]` using `[Cx,Cy]=chirp_pulse(parameters.pulsenpoints,parameters.Te, parameters.BW,parameters.nWURST,parameters.chirptype)`.
- Lines 68: computes `time_grid` using `time_grid=parameters.Te*ones(1,parameters.pulsenpoints)/parameters.pulsenpoints`.
- Lines 71: computes `gradient_amplitudes` using `gradient_amplitudes=parameters.Ge(1)*ones(parameters.pulsenpoints,1)`.
- Lines 74: computes `rho` using `rho=step(spin_system,Lx,parameters.rho0,pi/2)`.
- Lines 114: computes `Taq` using `Taq=parameters.npoints*parameters.deltat`.
- Lines 118: computes `fid` using `fid=zeros([parameters.npoints parameters.nloops],'like',1i)`.
- Lines 127: computes `EGp` using `EGp=L+parameters.Ga*G{1}`.
- Lines 128: computes `EGm` using `EGm=L-parameters.Ga*G{1}`.
- Lines 143: computes `fid(k,m)` using `fid(k,m)=parameters.coil'*rho`.

### Local helper functions

- Line 160: `grumble()` — `function grumble(spin_system,parameters,H,R,K,G,F)`.
  - Representative operation: `if ~ismember(spin_system.bas.formalism,{'sphten-liouv'})`.
  - Representative operation: `error('this function is only available for sphten-liouv formalism.')`.

## Outputs

- fid -free induction decay of the ultrafast NMR spectrum.

## Implementation structure

- Ultrafast multiple-quantum NMR, a literal implementation of Figure 1A
- from (http://dx.doi.org/10.1002/cphc.201800667). Syntax:
- fid=ufmq_nmr(spin_system,parameters,H,R,K,G,F)
- This sequence must be called from the imaging() context, which
- would provide H, R, K, G, and F. Parameters:
- parameters.spins nuclei on which the sequence runs
- parameters.dims size of the sample, m
- parameters.npts number of grid points
- parameters.npoints number of acquired points for each
- gradient readout
- parameters.nloops number of loop, where each loop consists of
- a positive and a negative readout

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `grumble()`, `operator()`, `polyadic()`, `speye()`, `ismember()`, `inflate()`, `chirp_pulse()`, `step()`, `evolution()`, `coherence()`, `shaped_pulse_xy()`, `gpuArray()`, `report()`, `num2str()`, `fid()`, `gather()`.
