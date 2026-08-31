# experiments/nmr_protein/hncoca.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/experiments/nmr_protein/hncoca.m`
- Signature: `fid=hncoca(spin_system,parameters,H,R,K)`
- Total lines: 272

## Purpose

Phase-sensitive HN(CO)CA pulse sequence from using the bidirectional propagation method described in The sequence is hard-wired to work on 1H,13C,15N proteins and uses PDB labels to select spins that will be affected by otherwise ideal pulses. F1 is N, F2 is CA, F3 is H. Syntax: fid=hncoca(spin_system,parameters,H,R,K)

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

- Lines 52-53: Consistency check; implemented by `grumble(spin_system,parameters,H,R,K)`.
- Lines 55-56: Compose Liouvillian; implemented by `L=H+1i*R+1i*K`.
- Lines 58-59: Coherent evolution timesteps; implemented by `t1.nsteps=parameters.npoints(1); t1.timestep=1/parameters.sweep(1)`.
- Lines 63-64: Initial condition -NH protons; implemented by `if ~isfield(parameters,'rho0')`.
- Lines 69-70: Detection state -NH protons; implemented by `if ~isfield(parameters,'coil')`.
- Lines 75-76: Spin indices; implemented by `COs=strcmp('C',spin_system.comp.labels)`.
- Lines 80-81: Pulse operators on 13CO carbons; implemented by `COp=operator(spin_system,'L+',find(COs))`.
- Lines 83-84: Pulse operators on 13CA carbons; implemented by `CAp=operator(spin_system,'L+',find(CAs))`.
- Lines 86-87: Pulse operators on all protons; implemented by `Hp=operator(spin_system,'L+','1H')`.
- Lines 89-90: Pulse operators on all nitrogens; implemented by `Np=operator(spin_system,'L+','15N')`.
- Lines 92-93: Cartesian pulse operators; implemented by `Hx=(Hp+Hp')/2; Hy=(Hp-Hp')/2i`.
- Lines 97-100: % Run the first half forward; implemented by `rho=step(spin_system,Hx,parameters.rho0,pi/2)`.
- Lines 99-100: Pulse on 1H; implemented by `rho=step(spin_system,Hx,parameters.rho0,pi/2)`.
- Lines 102-103: tau1 evolution; implemented by `rho=evolution(spin_system,L,[],rho,parameters.tau(1),1,'final')`.
- Lines 105-106: Inversion pulse on 1H, y inversion pulse on 15N; implemented by `rho=step(spin_system,Hx+Ny,rho,pi)`.
- Lines 111-112: Pulse on 15N, y pulse on 1H; implemented by `rho=step(spin_system,Hy+Nx,rho,pi/2)`.
- Lines 114-115: tau2 evolution; implemented by `rho=evolution(spin_system,L,[],rho,parameters.tau(2),1,'final')`.
- Lines 117-118: Coherence selection for States quadrature in F1; implemented by `rho_pos=coherence(spin_system,rho,{{find(Ns),+1}})`.

### Control flow inferred from the code

- Line 64: conditional branch on `~isfield(parameters,'rho0')`.
- Line 70: conditional branch on `~isfield(parameters,'coil')`.
- Line 223: `for` loop over `name_idx=1:numel(fid_names)`.

### Key state/data transformations

- Lines 56: computes `L` using `L=H+1i*R+1i*K`.
- Lines 59: computes `t1.nsteps` using `t1.nsteps=parameters.npoints(1); t1.timestep=1/parameters.sweep(1)`.
- Lines 60: computes `t2.nsteps` using `t2.nsteps=parameters.npoints(2); t2.timestep=1/parameters.sweep(2)`.
- Lines 61: computes `t3.nsteps` using `t3.nsteps=parameters.npoints(3); t3.timestep=1/parameters.sweep(3)`.
- Lines 65: computes `HNs` using `HNs=ismember(spin_system.comp.labels,{'H'})`.
- Lines 66: computes `parameters.rho0` using `parameters.rho0=state(spin_system,'Lz',find(HNs),'cheap')`.
- Lines 72: computes `parameters.coil` using `parameters.coil=state(spin_system,'L+',find(HNs),'cheap')`.
- Lines 76: computes `COs` using `COs=strcmp('C',spin_system.comp.labels)`.
- Lines 77: computes `CAs` using `CAs=strcmp('CA',spin_system.comp.labels)`.
- Lines 78: computes `Ns` using `Ns=strcmp('N',spin_system.comp.labels)`.
- Lines 81: computes `COp` using `COp=operator(spin_system,'L+',find(COs))`.
- Lines 84: computes `CAp` using `CAp=operator(spin_system,'L+',find(CAs))`.
- Lines 87: computes `Hp` using `Hp=operator(spin_system,'L+','1H')`.
- Lines 90: computes `Np` using `Np=operator(spin_system,'L+','15N')`.
- Lines 93: computes `Hx` using `Hx=(Hp+Hp')/2; Hy=(Hp-Hp')/2i`.
- Lines 94: computes `Nx` using `Nx=(Np+Np')/2; Ny=(Np-Np')/2i`.
- Lines 95: computes `COx` using `COx=(COp+COp')/2; CAx=(CAp+CAp')/2`.
- Lines 100: computes `rho` using `rho=step(spin_system,Hx,parameters.rho0,pi/2)`.

### Local helper functions

- Line 230: `grumble()` — `function grumble(spin_system,parameters,H,R,K)`.
  - Representative operation: `if ~ismember(spin_system.bas.formalism,{'sphten-liouv'})`.
  - Representative operation: `error('this function is only available for sphten-liouv formalism.')`.

## Parameters / inputs

- parameters.npoints -a vector of three integers giving the
- number of points in the three temporal
- dimensions, ordered as [t1 t2 t3].
- parameters.sweep -a vector of three real numbers giving
- the sweep widths in the three frequen-
- cy dimensions, ordered as [f1 f2 f3].
- parameters.tau -the four delays required for the ope-
- ration of the sequence (see the paper)
- in seconds. Reasonable values are
- [2.25e-3, 2.75e-3, 8.00e-3, 7.00e-3]
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

- Phase-sensitive HN(CO)CA pulse sequence from
- using the bidirectional propagation method described in
- The sequence is hard-wired to work on 1H,13C,15N proteins and uses
- PDB labels to select spins that will be affected by otherwise ideal
- pulses. F1 is N, F2 is CA, F3 is H. Syntax:
- fid=hncoca(spin_system,parameters,H,R,K)
- parameters.npoints -a vector of three integers giving the
- number of points in the three temporal
- dimensions, ordered as [t1 t2 t3].
- parameters.sweep -a vector of three real numbers giving
- the sweep widths in the three frequen-
- cy dimensions, ordered as [f1 f2 f3].

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `grumble()`, `isfield()`, `ismember()`, `state()`, `strcmp()`, `operator()`, `step()`, `evolution()`, `coherence()`, `decouple()`, `report()`, `stitch()`, `fieldnames()`, `ismatrix()`, `all()`, `isvector()`.
