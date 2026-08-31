# experiments/spen/spendosycosy.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/experiments/spen/spendosycosy.m`
- Signature: `fid=spendosycosy(spin_system,parameters,H,R,K,G,F)`
- Total lines: 321

## Purpose

Ultrafast 3D DOSY-COSY pulse sequence. Syntax: fid=spendosycosy(spin_system,parameters,H,R,K,G,F)

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
- Lines 80-81: Timestep for COSY; implemented by `timestep=1/parameters.sweep`.
- Lines 83-85: Get the chirp pulse; implemented by `[Cx,Cy]=chirp_pulse(parameters.pulsenpoints,parameters.Te,parameters.BW, parameters.smfactor,parameters.chirptype)`.
- Lines 87-88: Get the gradient amplitude list; implemented by `gradient_amplitudes=parameters.Ge*ones(1,parameters.pulsenpoints)`.
- Lines 90-91: Apply the first pulse; implemented by `rho=step(spin_system,Lx,parameters.rho0,pi/2)`.
- Lines 93-94: Select "+1" coherence; implemented by `rho=coherence(spin_system,rho,{{parameters.spins{1},+1}})`.
- Lines 96-97: Apply chirp pulse with a positive gradient; implemented by `time_grid=parameters.Te*ones(1,parameters.pulsenpoints)/parameters.pulsenpoints`.
- Lines 101-102: Select "-1" coherence; implemented by `rho=coherence(spin_system,rho,{{parameters.spins{1},-1}})`.
- Lines 104-105: Evolve under a positive gradient; implemented by `rho=step(spin_system,L+parameters.Ge*G{1},rho,parameters.Tau)`.
- Lines 107-108: Apply the second pulse; implemented by `rho=step(spin_system,Lx,rho,pi/2)`.
- Lines 110-111: Select "0" coherence; implemented by `rho=coherence(spin_system,rho,{{parameters.spins{1},0}})`.
- Lines 113-114: Run the diffusion time evolution; implemented by `rho=step(spin_system,L,rho,parameters.td-parameters.Tau-parameters.Te)`.
- Lines 116-117: Apply the third pulse; implemented by `rho=step(spin_system,Lx,rho,pi/2)`.
- Lines 122-124: Apply chirp pulse with a positive gradient; implemented by `rho=shaped_pulse_xy(spin_system,L,{Lx,Ly,G{1}},{Cx,Cy,+gradient_amplitudes}, time_grid,rho,'expv-pwc')`.
- Lines 132-133: Generate COSY stack; implemented by `cosy_stack=evolution(spin_system,L,[],rho,timestep,parameters.npoints2-1,'trajectory')`.
- Lines 135-136: Apply COSY pulse; implemented by `cosy_stack=step(spin_system,Lx,cosy_stack,pi/2)`.

### Control flow inferred from the code

- Line 154: conditional branch on `ismember('gpu',spin_system.sys.enable)`.
- Line 161: `for` loop over `m=1:parameters.nloops`.
- Line 173: `parfor` loop over `m=1:parameters.nloops`.
- Line 179: conditional branch on `ismember('gpu',spin_system.sys.enable)`.
- Line 189: `for` loop over `n=1:parameters.npoints1`.

### Key state/data transformations

- Lines 73: computes `L` using `L=H+F+1i*R+1i*K`.
- Lines 76: computes `Lp` using `Lp=operator(spin_system,'L+',parameters.spins{1})`.
- Lines 77: computes `Lx` using `Lx=kron(speye(prod(parameters.npts)),(Lp+Lp')/2)`.
- Lines 78: computes `Ly` using `Ly=kron(speye(prod(parameters.npts)),(Lp-Lp')/2i)`.
- Lines 81: computes `timestep` using `timestep=1/parameters.sweep`.
- Lines 84-85: computes `[Cx,Cy]` using `[Cx,Cy]=chirp_pulse(parameters.pulsenpoints,parameters.Te,parameters.BW, parameters.smfactor,parameters.chirptype)`.
- Lines 88: computes `gradient_amplitudes` using `gradient_amplitudes=parameters.Ge*ones(1,parameters.pulsenpoints)`.
- Lines 91: computes `rho` using `rho=step(spin_system,Lx,parameters.rho0,pi/2)`.
- Lines 97: computes `time_grid` using `time_grid=parameters.Te*ones(1,parameters.pulsenpoints)/parameters.pulsenpoints`.
- Lines 133: computes `cosy_stack` using `cosy_stack=evolution(spin_system,L,[],rho,timestep,parameters.npoints2-1,'trajectory')`.
- Lines 142: computes `Taq` using `Taq=parameters.npoints1*parameters.deltat`.
- Lines 146: computes `P1L` using `P1L=propagator(spin_system,L+parameters.Ga*G{1},parameters.npoints1*parameters.deltat)`.
- Lines 147: computes `P2L` using `P2L=propagator(spin_system,L-parameters.Ga*G{1},parameters.npoints1*parameters.deltat)`.
- Lines 148: computes `PL` using `PL=clean_up(spin_system,P2L*P1L,spin_system.tols.prop_chop); clear('P1L','P2L')`.
- Lines 151: computes `P` using `P=propagator(spin_system,L+parameters.Ga*G{1},parameters.deltat)`.
- Lines 160: computes `loop_stack` using `loop_stack=cell(1,parameters.nloops)`.
- Lines 162: computes `loop_stack{m}` using `loop_stack{m}=gather(cosy_stack)`.
- Lines 167-169: computes `fid` using `fid=zeros([parameters.npoints1 parameters.npoints2 parameters.nloops],'like',1i)`.

### Local helper functions

- Line 202: `grumble()` — `function grumble(spin_system,parameters,H,R,K,G,F)`.
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
- parameters.smfactor smoothing factor for the pulse
- parameters.Te duration of the pulse
- parameters.Tau duration extra dephasing gradient
- parameters.td diffusion delay, at least
- parameters.Tau+parameters.Te
- parameters.BW bandwidth of the pulse
- parameters.Ge encoding gradient in T/m
- parameters.Gp coherence selection gradient in T/m
- parameters.Tp duration of the coherence selection gradient
- parameters.chirptype can be 'wurst' or 'smoothed'
- H Fokker-Planck Hamiltonian
- R Fokker-Planck relaxation superoperator
- K Fokker-Planck kinetics superoperator
- G Fokker-Planck gradient superoperators
- F Fokker-Planck diffusion and flow superoperator

## Outputs

- fid -SPENDOSYCOSY free induction decay as a 3D array
- Note: the last five parameters are built automatically by the imaging
- context function.

## Implementation structure

- Ultrafast 3D DOSY-COSY pulse sequence. Syntax:
- fid=spendosycosy(spin_system,parameters,H,R,K,G,F)
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

- Called routines detected from the main body: `grumble()`, `operator()`, `speye()`, `chirp_pulse()`, `step()`, `coherence()`, `shaped_pulse_xy()`, `evolution()`, `propagator()`, `clean_up()`, `clear()`, `ismember()`, `gpuArray()`, `report()`, `gather()`, `local_fid()`.
