# experiments/nmr_protein/hnco.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/experiments/nmr_protein/hnco.m`
- Signature: `fid=hnco(spin_system,parameters,H,R,K)`
- Total lines: 254

## Purpose

Phase-sensitive HNCO pulse sequence from using the bidirectional propagation method described in The sequence is hard-wired to work on 1H,13C,15N proteins and uses PDB labels to select spins that will be affected by otherwise ideal pulses. F1 is N, F2 is CO, F3 is H. Syntax: fid=hnco(spin_system,parameters,H,R,K)

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

- Lines 55-56: Consistency check; implemented by `grumble(spin_system,parameters,H,R,K)`.
- Lines 58-59: Compose Liouvillian; implemented by `L=H+1i*R+1i*K`.
- Lines 61-62: Coherent evolution timesteps; implemented by `t1.nsteps=parameters.npoints(1); t1.timestep=1./parameters.sweep(1)`.
- Lines 66-67: Spin indices; implemented by `HNs=strcmp('H',spin_system.comp.labels)`.
- Lines 72-73: Initial condition -NH protons; implemented by `rho=state(spin_system,'Lz',find(HNs),'cheap')`.
- Lines 75-76: Detection state -NH protons; implemented by `coil=state(spin_system,'L+',find(HNs),'cheap')`.
- Lines 78-79: Pulse operators on 13CO carbons; implemented by `COp=operator(spin_system,'L+',find(COs))`.
- Lines 81-82: Pulse operators on 13CA carbons; implemented by `CAp=operator(spin_system,'L+',find(CAs))`.
- Lines 84-85: Pulse operators on all protons; implemented by `Hp=operator(spin_system,'L+','1H')`.
- Lines 87-88: Pulse operators on all nitrogens; implemented by `Np=operator(spin_system,'L+','15N')`.
- Lines 90-91: Cartesian pulse operators; implemented by `Nx=(Np+Np')/2; COx=(COp+COp')/2`.
- Lines 95-98: % Run the first half forward; implemented by `rho=step(spin_system,Hx,rho,pi/2)`.
- Lines 97-98: Pulse on 1H; implemented by `rho=step(spin_system,Hx,rho,pi/2)`.
- Lines 100-101: tau1 evolution; implemented by `rho=evolution(spin_system,L,[],rho,parameters.tau(1),1,'final')`.
- Lines 103-104: Inversion pulses on 1H and 15N; implemented by `rho=step(spin_system,Hx+Nx,rho,pi)`.
- Lines 109-110: Pulse on 15N, y pulse on 1H; implemented by `rho=step(spin_system,Hy+Nx,rho,pi/2)`.
- Lines 112-113: Coherence selection for States quadrature in F1; implemented by `rho_pos=coherence(spin_system,rho,{{find(Ns),+1}})`.
- Lines 116-117: t1 evolution; implemented by `rho_stack_pos=evolution(spin_system,L,[],rho_pos,t1.timestep/2,t1.nsteps-1,'trajectory')`.

### Control flow inferred from the code

- Line 121: conditional branch on `parameters.f1_decouple`.
- Line 200: `for` loop over `name_idx=1:numel(fid_names)`.

### Key state/data transformations

- Lines 59: computes `L` using `L=H+1i*R+1i*K`.
- Lines 62: computes `t1.nsteps` using `t1.nsteps=parameters.npoints(1); t1.timestep=1./parameters.sweep(1)`.
- Lines 63: computes `t2.nsteps` using `t2.nsteps=parameters.npoints(2); t2.timestep=1./parameters.sweep(2)`.
- Lines 64: computes `t3.nsteps` using `t3.nsteps=parameters.npoints(3); t3.timestep=1./parameters.sweep(3)`.
- Lines 67: computes `HNs` using `HNs=strcmp('H',spin_system.comp.labels)`.
- Lines 68: computes `COs` using `COs=strcmp('C',spin_system.comp.labels)`.
- Lines 69: computes `CAs` using `CAs=strcmp('CA',spin_system.comp.labels)`.
- Lines 70: computes `Ns` using `Ns=strcmp('N',spin_system.comp.labels)`.
- Lines 73: computes `rho` using `rho=state(spin_system,'Lz',find(HNs),'cheap')`.
- Lines 76: computes `coil` using `coil=state(spin_system,'L+',find(HNs),'cheap')`.
- Lines 79: computes `COp` using `COp=operator(spin_system,'L+',find(COs))`.
- Lines 82: computes `CAp` using `CAp=operator(spin_system,'L+',find(CAs))`.
- Lines 85: computes `Hp` using `Hp=operator(spin_system,'L+','1H')`.
- Lines 88: computes `Np` using `Np=operator(spin_system,'L+','15N')`.
- Lines 91: computes `Nx` using `Nx=(Np+Np')/2; COx=(COp+COp')/2`.
- Lines 92: computes `CAx` using `CAx=(CAp+CAp')/2; Hx=(Hp+Hp')/2`.
- Lines 93: computes `Hy` using `Hy=(Hp-Hp')/2i`.
- Lines 113: computes `rho_pos` using `rho_pos=coherence(spin_system,rho,{{find(Ns),+1}})`.

### Local helper functions

- Line 207: `grumble()` — `function grumble(spin_system,parameters,H,R,K)`.
  - Representative operation: `if ~ismember(spin_system.bas.formalism,{'sphten-liouv'})`.
  - Representative operation: `error('this function is only available for sphten-liouv formalisms.')`.

## Parameters / inputs

- parameters.npoints -a vector of three integers giving the
- number of points in the three temporal
- dimensions, ordered as [t1 t2 t3].
- parameters.sweep -a vector of three real numbers giving
- the sweep widths in the three frequen-
- cy dimensions, ordered as [f1 f2 f3].
- parameters.tau -the three delays required for the ope-
- ration of the sequence (see the paper)
- in seconds. Reasonable values are
- [2.25e-3, 14e-3, 4e-3]
- parameters.f1_decouple -logical switch controlling proton de-
- coupling during the T1 period.
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

- Phase-sensitive HNCO pulse sequence from
- using the bidirectional propagation method described in
- The sequence is hard-wired to work on 1H,13C,15N proteins and uses
- PDB labels to select spins that will be affected by otherwise ideal
- pulses. F1 is N, F2 is CO, F3 is H. Syntax:
- fid=hnco(spin_system,parameters,H,R,K)
- parameters.npoints -a vector of three integers giving the
- number of points in the three temporal
- dimensions, ordered as [t1 t2 t3].
- parameters.sweep -a vector of three real numbers giving
- the sweep widths in the three frequen-
- cy dimensions, ordered as [f1 f2 f3].

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `grumble()`, `strcmp()`, `state()`, `operator()`, `step()`, `evolution()`, `coherence()`, `decouple()`, `report()`, `stitch()`, `fieldnames()`, `ismember()`, `ismatrix()`, `all()`, `isfield()`, `isvector()`.
