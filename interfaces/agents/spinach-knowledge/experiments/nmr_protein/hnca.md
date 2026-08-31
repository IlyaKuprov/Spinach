# experiments/nmr_protein/hnca.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/experiments/nmr_protein/hnca.m`
- Signature: `fid=hnca(spin_system,parameters,H,R,K)`
- Total lines: 221

## Purpose

Protein-specific HNCA experiment (Figure 7.31a of "Protein NMR Spectroscopy", 2nd edition) using pre-set values of J-couplings used in the magnetisation transfer stages. The simulation uses the bidirectional propagation method described in The sequence is hard-wired to work on 1H,13C,15N proteins and uses PDB labels to select spins that will be affected by otherwise ideal pulses. F1 is 15N, F2 is 13C, F3 is 1H. Synta

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

- Lines 46-47: Check consistency; implemented by `grumble(spin_system,parameters,H,R,K)`.
- Lines 49-50: Compose Liouvillian; implemented by `L=H+1i*R+1i*K`.
- Lines 52-53: Coherent evolution timesteps; implemented by `t1.nsteps=parameters.npoints(1); t1.timestep=1./parameters.sweep(1)`.
- Lines 57-58: Hard-coded J-couplings and delays; implemented by `J_nh=92; tau=abs(1/(4*J_nh))`.
- Lines 61-62: Initial condition -NH protons; implemented by `if ~isfield(parameters,'rho0')`.
- Lines 67-68: Detection state -NH protons; implemented by `if ~isfield(parameters,'coil')`.
- Lines 73-74: Pulse operators all protons; implemented by `Hp=operator(spin_system,'L+','1H')`.
- Lines 77-78: Pulse operator on CA carbons; implemented by `CAs=strcmp('CA',spin_system.comp.labels)`.
- Lines 82-83: Pulse operators all nitrogens; implemented by `Np=operator(spin_system,'L+','15N')`.
- Lines 86-87: Pulse operator on CO carbons; implemented by `COs=strcmp('C',spin_system.comp.labels)`.
- Lines 91-94: % Run the first half forward; implemented by `rho=step(spin_system,Hx,parameters.rho0,pi/2)`.
- Lines 93-94: Pulse on 1H; implemented by `rho=step(spin_system,Hx,parameters.rho0,pi/2)`.
- Lines 96-97: Coherence selection; implemented by `rho=coherence(spin_system,rho,{{'1H',1}})`.
- Lines 99-100: tau evolution; implemented by `rho=evolution(spin_system,L,[],rho,tau,1,'final')`.
- Lines 102-103: Inversion pulses on 1H and 15N; implemented by `rho=step(spin_system,Hx+Nx,rho,pi)`.
- Lines 108-109: Pulse on 15N, y pulse on 1H; implemented by `rho=step(spin_system,Hy+Nx,rho,pi/2)`.
- Lines 111-112: Coherence selection for States quadrature in F1; implemented by `rho_pos=coherence(spin_system,rho,{{'15N',+1}})`.
- Lines 115-116: t1 evolution; implemented by `rho_stack_pos=evolution(spin_system,L,[],rho_pos,t1.timestep/2,t1.nsteps-1,'trajectory')`.

### Control flow inferred from the code

- Line 62: conditional branch on `~isfield(parameters,'rho0')`.
- Line 68: conditional branch on `~isfield(parameters,'coil')`.
- Line 179: `for` loop over `name_idx=1:numel(fid_names)`.

### Key state/data transformations

- Lines 50: computes `L` using `L=H+1i*R+1i*K`.
- Lines 53: computes `t1.nsteps` using `t1.nsteps=parameters.npoints(1); t1.timestep=1./parameters.sweep(1)`.
- Lines 54: computes `t2.nsteps` using `t2.nsteps=parameters.npoints(2); t2.timestep=1./parameters.sweep(2)`.
- Lines 55: computes `t3.nsteps` using `t3.nsteps=parameters.npoints(3); t3.timestep=1./parameters.sweep(3)`.
- Lines 58: computes `J_nh` using `J_nh=92; tau=abs(1/(4*J_nh))`.
- Lines 59: computes `J_nca` using `J_nca=11.5; delta=abs(1/(4*J_nca))`.
- Lines 63: computes `HNs` using `HNs=ismember(spin_system.comp.labels,{'H'})`.
- Lines 64: computes `parameters.rho0` using `parameters.rho0=state(spin_system,'Lz',find(HNs),'cheap')`.
- Lines 70: computes `parameters.coil` using `parameters.coil=state(spin_system,'L+',find(HNs),'cheap')`.
- Lines 74: computes `Hp` using `Hp=operator(spin_system,'L+','1H')`.
- Lines 75: computes `Hx` using `Hx=(Hp+Hp')/2; Hy=(Hp-Hp')/2i`.
- Lines 78: computes `CAs` using `CAs=strcmp('CA',spin_system.comp.labels)`.
- Lines 79: computes `CAp` using `CAp=operator(spin_system,'L+',find(CAs))`.
- Lines 80: computes `CAx` using `CAx=(CAp+CAp')/2`.
- Lines 83: computes `Np` using `Np=operator(spin_system,'L+','15N')`.
- Lines 84: computes `Nx` using `Nx=(Np+Np')/2`.
- Lines 87: computes `COs` using `COs=strcmp('C',spin_system.comp.labels)`.
- Lines 88: computes `COp` using `COp=operator(spin_system,'L+',find(COs))`.

### Local helper functions

- Line 186: `grumble()` — `function grumble(spin_system,parameters,H,R,K)`.
  - Representative operation: `if ~ismember(spin_system.bas.formalism,{'sphten-liouv'})`.
  - Representative operation: `error('this function is only available for sphten-liouv formalisms.')`.

## Parameters / inputs

- parameters.npoints -a vector of three integers giving the
- number of points in the three temporal
- dimensions, ordered as [t1 t2 t3].
- parameters.sweep -a vector of three real numbers giving
- the sweep widths in the three frequen-
- cy dimensions, ordered as [f1 f2 f3].
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

- Protein-specific HNCA experiment (Figure 7.31a of "Protein NMR
- Spectroscopy", 2nd edition) using pre-set values of J-couplings
- used in the magnetisation transfer stages. The simulation uses
- the bidirectional propagation method described in
- The sequence is hard-wired to work on 1H,13C,15N proteins and uses
- PDB labels to select spins that will be affected by otherwise ideal
- pulses. F1 is 15N, F2 is 13C, F3 is 1H. Syntax:
- fid=hnca(spin_system,parameters,H,R,K)
- parameters.npoints -a vector of three integers giving the
- number of points in the three temporal
- dimensions, ordered as [t1 t2 t3].
- parameters.sweep -a vector of three real numbers giving

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `grumble()`, `isfield()`, `ismember()`, `state()`, `operator()`, `strcmp()`, `step()`, `coherence()`, `evolution()`, `decouple()`, `stitch()`, `fieldnames()`, `ismatrix()`, `all()`, `isvector()`, `any()`.
