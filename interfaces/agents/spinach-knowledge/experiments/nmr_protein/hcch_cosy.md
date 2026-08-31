# experiments/nmr_protein/hcch_cosy.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/experiments/nmr_protein/hcch_cosy.m`
- Signature: `fid=hcch_cosy(spin_system,parameters,H,R,K)`
- Total lines: 287

## Purpose

HCCH-COSY pulse sequence from Figure 7.26a of Protein NMR Spectroscopy (2nd edition) using the bidirectional propagation method described in The sequence is hard-wired to work on 1H,13C proteins and uses PDB la- bels to select spins that will be affected by otherwise ideal pulses. F1 is 1H, F2 is C, F3 is H. Syntax: fid=hcch_cosy(spin_system,parameters,H,R,K)

## Physical / mathematical content

- Protein triple-resonance sequence implementations. They orchestrate heteronuclear coherence transfers across biomolecular spin networks while preserving phase and acquisition conventions.
- Orientation or trajectory averaging is performed numerically, so grid design, weights, and integration error control matter directly to accuracy and runtime.

## Numerical / algorithmic content

- Time propagation is explicit. In Spinach this usually means repeated application of matrix exponentials or propagator factorizations to density operators or state vectors in Hilbert/Liouville/Fokker-Planck space.
- Numerical integration over angles or geometry is part of the implementation, so point placement and weights are as important as the local Hamiltonian calculations.
- The file contains an explicit `grumble(...)` validator, which is Spinach convention for front-loading dimension, type, and regime checks before expensive linear-algebra work begins.
- The file also defines local helper function(s): `grumble()`. This usually means the public entry point is supported by tightly coupled validation or helper logic kept private to the file.

## Code-derived implementation details

### Comment-guided execution stages

- Lines 56-57: Consistency check; implemented by `grumble(spin_system,parameters,H,R,K)`.
- Lines 59-60: Compose Liouvillian; implemented by `L=H+1i*R+1i*K`.
- Lines 62-63: Coherent evolution timesteps; implemented by `t1.nsteps=parameters.npoints(1); t1.timestep=1./parameters.sweep(1)`.
- Lines 67-68: J-coupling evolution time; implemented by `tau_ch=abs(1/(4*parameters.J_ch))`.
- Lines 70-71: Coherence transfer delays; implemented by `tau_cc=abs(1/(8*parameters.J_cc))`.
- Lines 74-75: Initial condition; implemented by `rho0=state(spin_system,'Lz','1H','cheap')`.
- Lines 77-78: Detection state; implemented by `coil=state(spin_system,'L+','1H','cheap')`.
- Lines 80-81: Pulse operators all protons; implemented by `Hp=operator(spin_system,'L+','1H'); Hx=(Hp+Hp')/2`.
- Lines 83-84: Pulse operators all carbons; implemented by `Cp=operator(spin_system,'L+','13C'); Cx=(Cp+Cp')/2`.
- Lines 86-87: Selective pulse on CO carbons; implemented by `COs=strcmp('C',spin_system.comp.labels)`.
- Lines 91-94: % Forward sim from rho0 up to t2 period; implemented by `rho=step(spin_system,Hx,rho0,pi/2)`.
- Lines 93-94: Pulse on 1H; implemented by `rho=step(spin_system,Hx,rho0,pi/2)`.
- Lines 96-97: Coherence selection for States quadrature in F1; implemented by `rho_pos=coherence(spin_system,rho,{{'1H',+1}})`.
- Lines 100-101: t1 evolution; implemented by `rho_stack_pos=evolution(spin_system,L,[],rho_pos,t1.timestep/2,t1.nsteps-1,'trajectory')`.
- Lines 104-105: Inversion pulse on 13C; implemented by `rho_stack_pos=step(spin_system,Cx,rho_stack_pos,pi)`.
- Lines 108-109: t1 rest of the evolution; implemented by `rho_stack_pos=evolution(spin_system,L,[],rho_stack_pos,t1.timestep/2,t1.nsteps-1,'refocus')`.
- Lines 112-113: tau evolution; implemented by `rho_stack_pos=evolution(spin_system,L,[],rho_stack_pos,tau_ch,1,'final')`.
- Lines 116-117: Inversion pulses on 1H and 13C; implemented by `rho_stack_pos=step(spin_system,Hx+Cx,rho_stack_pos,pi)`.

### Control flow inferred from the code

- Line 194: `for` loop over `name_idx=1:numel(fid_names)`.

### Key state/data transformations

- Lines 60: computes `L` using `L=H+1i*R+1i*K`.
- Lines 63: computes `t1.nsteps` using `t1.nsteps=parameters.npoints(1); t1.timestep=1./parameters.sweep(1)`.
- Lines 64: computes `t2.nsteps` using `t2.nsteps=parameters.npoints(2); t2.timestep=1./parameters.sweep(2)`.
- Lines 65: computes `t3.nsteps` using `t3.nsteps=parameters.npoints(3); t3.timestep=1./parameters.sweep(3)`.
- Lines 68: computes `tau_ch` using `tau_ch=abs(1/(4*parameters.J_ch))`.
- Lines 71: computes `tau_cc` using `tau_cc=abs(1/(8*parameters.J_cc))`.
- Lines 72: computes `DELTA` using `DELTA=tau_cc-parameters.delta`.
- Lines 75: computes `rho0` using `rho0=state(spin_system,'Lz','1H','cheap')`.
- Lines 78: computes `coil` using `coil=state(spin_system,'L+','1H','cheap')`.
- Lines 81: computes `Hp` using `Hp=operator(spin_system,'L+','1H'); Hx=(Hp+Hp')/2`.
- Lines 84: computes `Cp` using `Cp=operator(spin_system,'L+','13C'); Cx=(Cp+Cp')/2`.
- Lines 87: computes `COs` using `COs=strcmp('C',spin_system.comp.labels)`.
- Lines 88: computes `COp` using `COp=operator(spin_system,'L+',find(COs))`.
- Lines 89: computes `COx` using `COx=(COp+COp')/2`.
- Lines 94: computes `rho` using `rho=step(spin_system,Hx,rho0,pi/2)`.
- Lines 97: computes `rho_pos` using `rho_pos=coherence(spin_system,rho,{{'1H',+1}})`.
- Lines 98: computes `rho_neg` using `rho_neg=coherence(spin_system,rho,{{'1H',-1}})`.
- Lines 101: computes `rho_stack_pos` using `rho_stack_pos=evolution(spin_system,L,[],rho_pos,t1.timestep/2,t1.nsteps-1,'trajectory')`.

### Local helper functions

- Line 201: `grumble()` — `function grumble(spin_system,parameters,H,R,K)`.
  - Representative operation: `if ~ismember(spin_system.bas.formalism,{'sphten-liouv'})`.
  - Representative operation: `error('this function is only available for sphten-liouv formalism.')`.

## Parameters / inputs

- parameters.npoints -a vector of three integers giving the
- number of points in the three temporal
- dimensions, ordered as [t1 t2 t3].
- parameters.sweep -a vector of three real numbers giving
- the sweep widths in the three frequen-
- cy dimensions, ordered as [f1 f2 f3].
- parameters.J_cc -13C-13C J-coupling to be used for mag-
- netisation transfer, typically 35 Hz
- parameters.J_ch -1H-13C J-coupling to be used for mag-
- netisation transfer, typically 140 Hz
- parameters.delta -evolution delay, see the pulse sequence
- diagram, typically 1.1e-3 seconds.
- parameters.decouple_f3 -list of spins to be decoupled during
- the detection period, typically {'13C'}
- H -Hamiltonian matrix, received from context function
- R -relaxation superoperator, received from context function
- K -kinetics superoperator, received from context function

## Outputs

- fid -a structure with four fields: fid.pos_pos, fid.pos_neg,
- fid.neg_pos, fid.neg_neg that are used in the subsequ-
- ent States quadrature processing
- Note: spin labels must be set to PDB atom IDs ('CA', 'HA', etc.) in
- sys.labels for this sequence to work properly.

## Implementation structure

- HCCH-COSY pulse sequence from Figure 7.26a of Protein NMR Spectroscopy
- (2nd edition) using the bidirectional propagation method described in
- The sequence is hard-wired to work on 1H,13C proteins and uses PDB la-
- bels to select spins that will be affected by otherwise ideal pulses.
- F1 is 1H, F2 is C, F3 is H. Syntax:
- fid=hcch_cosy(spin_system,parameters,H,R,K)
- parameters.npoints -a vector of three integers giving the
- number of points in the three temporal
- dimensions, ordered as [t1 t2 t3].
- parameters.sweep -a vector of three real numbers giving
- the sweep widths in the three frequen-
- cy dimensions, ordered as [f1 f2 f3].

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `grumble()`, `state()`, `operator()`, `strcmp()`, `step()`, `coherence()`, `evolution()`, `decouple()`, `stitch()`, `fieldnames()`, `ismember()`, `ismatrix()`, `all()`, `isfield()`, `isvector()`, `any()`.
