# experiments/spen/spendosy.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/experiments/spen/spendosy.m`
- Signature: `fid=spendosy(spin_system,parameters,H,R,K,G,F)`
- Total lines: 285

## Purpose

Ultrafast DOSY pulse sequence. Syntax: fid=spendosy(spin_system,parameters,H,R,K,G,F)

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

- Lines 69-70: Check consistency; implemented by `grumble(spin_system,parameters,H,R,K,G,F)`.
- Lines 72-73: Compose Liouvillian; implemented by `L=H+F+1i*R+1i*K`.
- Lines 75-76: Get pulse operators; implemented by `Lp=operator(spin_system,'L+',parameters.spins{1})`.
- Lines 80-82: Get the chirp pulse; implemented by `[Cx,Cy]=chirp_pulse(parameters.pulsenpoints,parameters.Te,parameters.BW, parameters.smfactor,parameters.chirptype)`.
- Lines 84-85: Get the gradient amplitude list; implemented by `gradient_amplitudes=parameters.Ge*ones(parameters.pulsenpoints,1)`.
- Lines 87-88: Apply the first pulse; implemented by `rho=step(spin_system,Lx,parameters.rho0,pi/2)`.
- Lines 90-91: Select "+1" coherence; implemented by `rho=coherence(spin_system,rho,{{parameters.spins{1},+1}})`.
- Lines 93-94: Apply chirp pulse with a positive gradient; implemented by `time_grid=parameters.Te*ones(1,parameters.pulsenpoints)/parameters.pulsenpoints`.
- Lines 98-99: Select "-1" coherence; implemented by `rho=coherence(spin_system,rho,{{parameters.spins{1},-1}})`.
- Lines 101-102: Evolve under a positive gradient; implemented by `rho=step(spin_system,L+parameters.Ge*G{1},rho,parameters.Tau)`.
- Lines 104-105: Apply the second pulse; implemented by `rho=step(spin_system,Lx,rho,pi/2)`.
- Lines 107-108: Select "0" coherence; implemented by `rho=coherence(spin_system,rho,{{parameters.spins{1},0}})`.
- Lines 110-111: Run the diffusion time evolution; implemented by `rho=step(spin_system,L,rho,parameters.td-parameters.Tau-parameters.Te)`.
- Lines 113-114: Apply the third pulse; implemented by `rho=step(spin_system,Lx,rho,pi/2)`.
- Lines 119-121: Apply chirp pulse with a positive gradient; implemented by `rho=shaped_pulse_xy(spin_system,L,{Lx,Ly,G{1}},{Cx,Cy,+gradient_amplitudes}, time_grid,rho,'expv-pwc')`.
- Lines 129-130: Apply prephasing; implemented by `Taq=parameters.npoints*parameters.deltat`.
- Lines 133-134: Build whole-loop propagators; implemented by `P1L=propagator(spin_system,L+parameters.Ga*G{1},parameters.npoints*parameters.deltat)`.
- Lines 138-139: Build the intra-loop propagator; implemented by `P=propagator(spin_system,L+parameters.Ga*G{1},parameters.deltat)`.

### Control flow inferred from the code

- Line 142: conditional branch on `ismember('gpu',spin_system.sys.enable)`.
- Line 153: `for` loop over `m=1:parameters.nloops`.
- Line 163: `parfor` loop over `m=1:parameters.nloops`.
- Line 166: `for` loop over `n=1:parameters.npoints`.

### Key state/data transformations

- Lines 73: computes `L` using `L=H+F+1i*R+1i*K`.
- Lines 76: computes `Lp` using `Lp=operator(spin_system,'L+',parameters.spins{1})`.
- Lines 77: computes `Lx` using `Lx=kron(speye(prod(parameters.npts)),(Lp+Lp')/2)`.
- Lines 78: computes `Ly` using `Ly=kron(speye(prod(parameters.npts)),(Lp-Lp')/2i)`.
- Lines 81-82: computes `[Cx,Cy]` using `[Cx,Cy]=chirp_pulse(parameters.pulsenpoints,parameters.Te,parameters.BW, parameters.smfactor,parameters.chirptype)`.
- Lines 85: computes `gradient_amplitudes` using `gradient_amplitudes=parameters.Ge*ones(parameters.pulsenpoints,1)`.
- Lines 88: computes `rho` using `rho=step(spin_system,Lx,parameters.rho0,pi/2)`.
- Lines 94: computes `time_grid` using `time_grid=parameters.Te*ones(1,parameters.pulsenpoints)/parameters.pulsenpoints`.
- Lines 130: computes `Taq` using `Taq=parameters.npoints*parameters.deltat`.
- Lines 134: computes `P1L` using `P1L=propagator(spin_system,L+parameters.Ga*G{1},parameters.npoints*parameters.deltat)`.
- Lines 135: computes `P2L` using `P2L=propagator(spin_system,L-parameters.Ga*G{1},parameters.npoints*parameters.deltat)`.
- Lines 136: computes `PL` using `PL=clean_up(spin_system,P2L*P1L,spin_system.tols.prop_chop); clear('P1L','P2L')`.
- Lines 139: computes `P` using `P=propagator(spin_system,L+parameters.Ga*G{1},parameters.deltat)`.
- Lines 145: computes `coil` using `coil=gpuArray(parameters.coil)`.
- Lines 152: computes `rho_stack` using `rho_stack=cell(1,parameters.nloops)`.
- Lines 154: computes `rho_stack{m}` using `rho_stack{m}=gather(rho)`.
- Lines 159: computes `fid` using `fid=zeros([parameters.npoints parameters.nloops],'like',1i)`.
- Lines 165: computes `local_fid` using `local_fid=zeros(parameters.npoints,1)`.

### Local helper functions

- Line 176: `grumble()` — `function grumble(spin_system,parameters,H,R,K,G,F)`.
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
- parameters.offset offset
- parameters.cond bondary conditions
- parameters.Ga acquisition gradient in T/m
- parameters.pulsenpoints number of points in the pulse shape
- parameters.smfactor smoothing factor for the pulse
- parameters.Te duration of the pulse
- parameters.Tau duration extra dephasing gradient
- parameters.BW bandwidth of the pulse
- parameters.Ge encoding gradient in T/m
- parameters.chirptype can be 'wurst' or 'smoothed'
- parameters.td diffusion delay, at least
- parameters.Tau+parameters.Te
- H Fokker-Planck Hamiltonian
- R Fokker-Planck relaxation superoperator
- K Fokker-Planck kinetics superoperator
- G Fokker-Planck gradient superoperators
- F Fokker-Planck diffusion and flow superoperator

## Outputs

- fid -free induction decay
- Note: the last five parameters are built automatically by the imaging
- context function.

## Implementation structure

- Ultrafast DOSY pulse sequence. Syntax:
- fid=spendosy(spin_system,parameters,H,R,K,G,F)
- parameters.dims size of the sample in m
- parameters.npts number of spin packets
- parameters.spins nuclei on which the sequence runs
- parameters.deltat timestep for acquisition
- parameters.npoints number of acquired points for each
- gradient readout
- parameters.nloops number of loop, where each loop consists of
- a positive and a negative readout
- parameters.offset offset
- parameters.cond bondary conditions

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `grumble()`, `operator()`, `speye()`, `chirp_pulse()`, `step()`, `coherence()`, `shaped_pulse_xy()`, `propagator()`, `clean_up()`, `clear()`, `ismember()`, `gpuArray()`, `report()`, `gather()`, `local_fid()`, `fid()`.
