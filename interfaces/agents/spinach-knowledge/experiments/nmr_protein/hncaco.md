# experiments/nmr_protein/hncaco.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/experiments/nmr_protein/hncaco.m`
- Signature: `fid=hncaco(spin_system,parameters,H,R,K)`
- Total lines: 274

## Purpose

Protein-specific HN(CA)CO experiment (Figure 7.41 of "Protein NMR Spectroscopy", 2nd edition) using pre-set values of J-couplings used in the magnetisation transfer stages. The simulation uses the bidirectional propagation method described in The sequence is hard-wired to work on 1H,13C,15N proteins and uses PDB labels to select spins that will be affected by otherwise ideal pulses. F1 is 15N, F2 is 13C, F3 is 1H. Sy

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

- Lines 54-55: Check consistency; implemented by `grumble(spin_system,parameters,H,R,K)`.
- Lines 57-58: Compose Liouvillian; implemented by `L=H+1i*R+1i*K`.
- Lines 60-61: Coherent evolution timesteps; implemented by `t1.nsteps=parameters.npoints(1); t1.timestep=1./parameters.sweep(1)`.
- Lines 65-66: J-coupling evolution time; implemented by `tau=abs(1/(4*parameters.J_nh))`.
- Lines 69-70: Pulse operators all protons; implemented by `Hp=operator(spin_system,'L+','1H')`.
- Lines 73-74: Spin indices and pulse operator -only on CA carbons; implemented by `CAs=strcmp('CA',spin_system.comp.labels)`.
- Lines 78-79: Pulse operators all N; implemented by `Np=operator(spin_system,'L+','15N')`.
- Lines 82-83: Spin indices and pulse operator -only on CO carbons; implemented by `COs=strcmp('C',spin_system.comp.labels)`.
- Lines 87-88: Detection state -NH protons; implemented by `HNs=strcmp('H',spin_system.comp.labels)`.
- Lines 91-94: % Run the first half forward; implemented by `rho_pos=state(spin_system,'L+','15N')`.
- Lines 93-94: Start in N+/-to emulate the INEPT block; implemented by `rho_pos=state(spin_system,'L+','15N')`.
- Lines 97-98: Decouple protons; implemented by `L_decH=decouple(spin_system,L,[],{'1H'})`.
- Lines 100-101: First half of t1 evolution; implemented by `rho_stack_pos=evolution(spin_system,L_decH,[],rho_pos,t1.timestep/2,t1.nsteps-1,'trajectory')`.
- Lines 104-105: Inversion pulse on CO; implemented by `rho_stack_pos=step(spin_system,COx,rho_stack_pos,pi)`.
- Lines 108-109: T/2 evolution; implemented by `rho_stack_pos=evolution(spin_system,L_decH,[],rho_stack_pos,parameters.T/2,1,'final')`.
- Lines 112-113: Inversion pulses on N and CA; implemented by `rho_stack_pos=step(spin_system,Nx+CAx,rho_stack_pos,pi)`.
- Lines 120-121: -t1/2 evolution; implemented by `rho_stack_pos=evolution(spin_system,L_decH,[],rho_stack_pos,-t1.timestep/2,t1.nsteps-1,'refocus')`.
- Lines 124-125: Pulses on 15N and 13CA; implemented by `rho_stack_pos=step(spin_system,Nx+CAx,rho_stack_pos,pi/2)`.

### Control flow inferred from the code

- Line 206: `for` loop over `name_idx=1:numel(fid_names)`.

### Key state/data transformations

- Lines 58: computes `L` using `L=H+1i*R+1i*K`.
- Lines 61: computes `t1.nsteps` using `t1.nsteps=parameters.npoints(1); t1.timestep=1./parameters.sweep(1)`.
- Lines 62: computes `t2.nsteps` using `t2.nsteps=parameters.npoints(2); t2.timestep=1./parameters.sweep(2)`.
- Lines 63: computes `t3.nsteps` using `t3.nsteps=parameters.npoints(3); t3.timestep=1./parameters.sweep(3)`.
- Lines 66: computes `tau` using `tau=abs(1/(4*parameters.J_nh))`.
- Lines 67: computes `delta` using `delta=abs(1/(2*parameters.J_nh))`.
- Lines 70: computes `Hp` using `Hp=operator(spin_system,'L+','1H')`.
- Lines 71: computes `Hx` using `Hx=(Hp+Hp')/2; Hy=(Hp-Hp')/2i`.
- Lines 74: computes `CAs` using `CAs=strcmp('CA',spin_system.comp.labels)`.
- Lines 75: computes `CAp` using `CAp=operator(spin_system,'L+',find(CAs))`.
- Lines 76: computes `CAx` using `CAx=(CAp+CAp')/2; CAy=(CAp-CAp')/2i`.
- Lines 79: computes `Np` using `Np=operator(spin_system,'L+','15N')`.
- Lines 80: computes `Nx` using `Nx=(Np+Np')/2`.
- Lines 83: computes `COs` using `COs=strcmp('C',spin_system.comp.labels)`.
- Lines 84: computes `COp` using `COp=operator(spin_system,'L+',find(COs))`.
- Lines 85: computes `COx` using `COx=(COp+COp')/2`.
- Lines 88: computes `HNs` using `HNs=strcmp('H',spin_system.comp.labels)`.
- Lines 89: computes `coil` using `coil=state(spin_system,'L+',find(HNs))`.

### Local helper functions

- Line 213: `grumble()` — `function grumble(spin_system,parameters,H,R,K)`.
  - Representative operation: `if ~ismember(spin_system.bas.formalism,{'sphten-liouv'})`.
  - Representative operation: `error('this function is only available for sphten-liouv formalisms.')`.

## Parameters / inputs

- parameters.npoints -a vector of three integers giving the
- number of points in the three temporal
- dimensions, ordered as [t1 t2 t3].
- parameters.sweep -a vector of three real numbers giving
- the sweep widths in the three frequen-
- cy dimensions, ordered as [f1 f2 f3].
- parameters.J_nh -1H-15N J-coupling in Hz to be used for
- magnetisation transfer.
- parameters.T -evolution delay in the indirect 15N
- dimension, in seconds.
- parameters.delta2 -coherence transfer delay in seconds.
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

- Protein-specific HN(CA)CO experiment (Figure 7.41 of "Protein NMR
- Spectroscopy", 2nd edition) using pre-set values of J-couplings
- used in the magnetisation transfer stages. The simulation uses
- the bidirectional propagation method described in
- The sequence is hard-wired to work on 1H,13C,15N proteins and uses
- PDB labels to select spins that will be affected by otherwise ideal
- pulses. F1 is 15N, F2 is 13C, F3 is 1H. Syntax:
- fid=hncaco(spin_system,parameters,H,R,K)
- parameters.npoints -a vector of three integers giving the
- number of points in the three temporal
- dimensions, ordered as [t1 t2 t3].
- parameters.sweep -a vector of three real numbers giving

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `grumble()`, `operator()`, `strcmp()`, `state()`, `decouple()`, `evolution()`, `step()`, `coherence()`, `stitch()`, `fieldnames()`, `ismember()`, `ismatrix()`, `all()`, `isfield()`, `isvector()`, `any()`.
