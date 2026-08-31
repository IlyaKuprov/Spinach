# experiments/spen/spencosy.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/experiments/spen/spencosy.m`
- Signature: `fid=spencosy(spin_system,parameters,H,R,K,G,F)`
- Total lines: 245

## Purpose

Ultrafast COSY pulse sequence. Syntax: fid=spencosy(spin_system,parameters,H,R,K,G,F)

## Physical / mathematical content

- SPEN experiment implementations. These files combine shaped pulses, gradients, spatial encoding, and often diffusion-aware propagation.
- The file uses a Fokker-Planck-style enlarged state space in which spatial or orientational coordinates are promoted to extra dimensions and coupled to spin dynamics through differential operators.

## Numerical / algorithmic content

- Time propagation is explicit. In Spinach this usually means repeated application of matrix exponentials or propagator factorizations to density operators or state vectors in Hilbert/Liouville/Fokker-Planck space.
- The implementation explicitly addresses performance engineering through parallel or GPU execution, which matters because Spinach operators can become extremely large after basis expansion or powder/spatial lifting.
- The file contains an explicit `grumble(...)` validator, which is Spinach convention for front-loading dimension, type, and regime checks before expensive linear-algebra work begins.
- The file also defines local helper function(s): `grumble()`. This usually means the public entry point is supported by tightly coupled validation or helper logic kept private to the file.

## Code-derived implementation details

### Comment-guided execution stages

- Lines 64-65: Check consistency; implemented by `grumble(spin_system,parameters,H,R,K,G,F)`.
- Lines 67-68: Compose Liouvillian; implemented by `L=H+F+1i*R+1i*K`.
- Lines 70-71: Get pulse operators; implemented by `Lp=operator(spin_system,'L+',parameters.spins{1})`.
- Lines 75-77: Get chirp pulse waveform; implemented by `[Cx,Cy]=chirp_pulse(parameters.pulsenpoints,parameters.Te, parameters.BW,parameters.nWURST,'wurst')`.
- Lines 79-80: Get the gradient pulse waveform; implemented by `gradient_amplitudes=parameters.Ge*ones(parameters.pulsenpoints,1)`.
- Lines 82-83: Apply the first pulse; implemented by `rho=step(spin_system,Lx,parameters.rho0,pi/2)`.
- Lines 85-86: Apply the pair of chirp pulses with opposite gradients; implemented by `report(spin_system,'running chirps ')`.
- Lines 93-94: Apply the second pulse with coherence selection; implemented by `rho=step(spin_system,parameters.Gp*G{1},rho,parameters.Tp)`.
- Lines 98-99: Apply prephasing; implemented by `Taq=parameters.npoints*parameters.deltat`.
- Lines 102-103: Build whole-loop propagators; implemented by `P1L=propagator(spin_system,L+parameters.Ga*G{1},parameters.npoints*parameters.deltat)`.
- Lines 107-108: Build the intra-loop propagator; implemented by `P=propagator(spin_system,L+parameters.Ga*G{1},parameters.deltat)`.
- Lines 110-111: Move to the GPU if necessary; implemented by `if ismember('gpu',spin_system.sys.enable)`.
- Lines 119-120: Generate loop starts; implemented by `report(spin_system,'computing loop starts ')`.
- Lines 126-127: Preallocate the fid; implemented by `fid=zeros([parameters.npoints parameters.nloops],'like',1i)`.
- Lines 129-130: Propagate the system; implemented by `report(spin_system,'computing loop bodies ')`.

### Control flow inferred from the code

- Line 111: conditional branch on `ismember('gpu',spin_system.sys.enable)`.
- Line 122: `for` loop over `m=1:parameters.nloops`.
- Line 131: `parfor` loop over `m=1:parameters.nloops`.
- Line 134: `for` loop over `n=1:parameters.npoints`.

### Key state/data transformations

- Lines 68: computes `L` using `L=H+F+1i*R+1i*K`.
- Lines 71: computes `Lp` using `Lp=operator(spin_system,'L+',parameters.spins{1})`.
- Lines 72: computes `Lx` using `Lx=kron(speye(prod(parameters.npts)),(Lp+Lp')/2)`.
- Lines 73: computes `Ly` using `Ly=kron(speye(prod(parameters.npts)),(Lp-Lp')/2i)`.
- Lines 76-77: computes `[Cx,Cy]` using `[Cx,Cy]=chirp_pulse(parameters.pulsenpoints,parameters.Te, parameters.BW,parameters.nWURST,'wurst')`.
- Lines 80: computes `gradient_amplitudes` using `gradient_amplitudes=parameters.Ge*ones(parameters.pulsenpoints,1)`.
- Lines 83: computes `rho` using `rho=step(spin_system,Lx,parameters.rho0,pi/2)`.
- Lines 87: computes `time_grid` using `time_grid=parameters.Te*ones(1,parameters.pulsenpoints)/parameters.pulsenpoints`.
- Lines 99: computes `Taq` using `Taq=parameters.npoints*parameters.deltat`.
- Lines 103: computes `P1L` using `P1L=propagator(spin_system,L+parameters.Ga*G{1},parameters.npoints*parameters.deltat)`.
- Lines 104: computes `P2L` using `P2L=propagator(spin_system,L-parameters.Ga*G{1},parameters.npoints*parameters.deltat)`.
- Lines 105: computes `PL` using `PL=clean_up(spin_system,P2L*P1L,spin_system.tols.prop_chop); clear('P1L','P2L')`.
- Lines 108: computes `P` using `P=propagator(spin_system,L+parameters.Ga*G{1},parameters.deltat)`.
- Lines 114: computes `coil` using `coil=gpuArray(parameters.coil)`.
- Lines 121: computes `rho_stack` using `rho_stack=cell(1,parameters.nloops)`.
- Lines 123: computes `rho_stack{m}` using `rho_stack{m}=rho; rho=PL*rho`.
- Lines 127: computes `fid` using `fid=zeros([parameters.npoints parameters.nloops],'like',1i)`.
- Lines 133: computes `local_fid` using `local_fid=zeros(parameters.npoints,1)`.

### Local helper functions

- Line 144: `grumble()` — `function grumble(spin_system,parameters,H,R,K,G,F)`.
  - Representative operation: `if ~ismember(spin_system.bas.formalism,{'sphten-liouv'})`.
  - Representative operation: `error('this function is only available for sphten-liouv formalism.')`.

## Parameters / inputs

- parameters.dims size of the sample in m
- parameters.npts number of spin packets
- parameters.spins nuclei on which the sequence runs
- parameters.deltat timestep for acquisition
- parameters.npoints number of acquired points for each
- gradient readout
- parameters.nloops number of loop, where each loop consists of
- a positive and a negative readout
- parameters.Ga acquisition gradient in T/m
- parameters.pulsenpoints number of points in the pulse shape
- parameters.nWURST smoothing factor for the pulse
- parameters.Te duration of the pulse
- parameters.BW bandwidth of the pulse
- parameters.Ge encoding gradient in T/m
- parameters.Gp coherence selection gradient in T/m
- parameters.Tp duration of the coherence selection gradient
- parameters.D diffusion constant, m^2/s
- H Fokker-Planck Hamiltonian
- R Fokker-Planck relaxation superoperator
- K Fokker-Planck kinetics superoperator
- G Fokker-Planck gradient superoperators
- F Fokker-Planck diffusion and flow superoperator

## Outputs

- fid UFCOSY free induction decay
- Note: the last five parameters are built automatically by the imaging
- context function.

## Implementation structure

- Ultrafast COSY pulse sequence. Syntax:
- fid=spencosy(spin_system,parameters,H,R,K,G,F)
- parameters.dims size of the sample in m
- parameters.npts number of spin packets
- parameters.spins nuclei on which the sequence runs
- parameters.deltat timestep for acquisition
- parameters.npoints number of acquired points for each
- gradient readout
- parameters.nloops number of loop, where each loop consists of
- a positive and a negative readout
- parameters.Ga acquisition gradient in T/m
- parameters.pulsenpoints number of points in the pulse shape

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `grumble()`, `operator()`, `speye()`, `chirp_pulse()`, `step()`, `report()`, `shaped_pulse_xy()`, `propagator()`, `clean_up()`, `clear()`, `ismember()`, `gpuArray()`, `local_fid()`, `gather()`, `fid()`, `ismatrix()`.
