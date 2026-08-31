# experiments/nmr_protein/hcanh.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/experiments/nmr_protein/hcanh.m`
- Signature: `fid=hcanh(spin_system,parameters,H,R,K)`
- Total lines: 272

## Purpose

Protein-specific H(CA)NH experiment (Figure 7.37 of "Protein NMR Spectroscopy", 2nd edition) using pre-set values of J-couplings used in the magnetisation transfer stages. The simulation uses the bidirectional propagation method described in The sequence is hard-wired to work on 1H,13C,15N proteins and uses PDB labels to select spins that will be affected by otherwise ideal pulses. F1 is 1H, F2 is 15N, F3 is 1H. Synt

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

- Lines 50-51: Consistency check; implemented by `grumble(spin_system,parameters,H,R,K)`.
- Lines 53-54: Compose Liouvillian; implemented by `L=H+1i*R+1i*K`.
- Lines 56-57: Coherent evolution timesteps; implemented by `t1.nsteps=parameters.npoints(1); t1.timestep=1./parameters.sweep(1)`.
- Lines 61-62: Hard-coded J-couplings (Hz) and delays (seconds); implemented by `J_ch=140; J_nh=92`.
- Lines 70-71: Initial condition -CA protons; implemented by `if ~isfield(parameters,'rho0')`.
- Lines 76-77: Detection state -NH protons; implemented by `if ~isfield(parameters,'coil')`.
- Lines 82-83: Pulse operators all protons; implemented by `Hp=operator(spin_system,'L+',parameters.spins{1})`.
- Lines 86-87: Spin indices and pulse operator -only on CA carbons; implemented by `CAs=strcmp('CA',spin_system.comp.labels)`.
- Lines 91-92: Pulse operators all N; implemented by `Np=operator(spin_system,'L+','15N')`.
- Lines 95-96: Spin indices on CO carbons; implemented by `COs=strcmp('C',spin_system.comp.labels)`.
- Lines 98-101: % Forward sim from rho0 up to t2 period; implemented by `rho=step(spin_system,Hx,parameters.rho0,pi/2)`.
- Lines 100-101: Pulse on 1H; implemented by `rho=step(spin_system,Hx,parameters.rho0,pi/2)`.
- Lines 103-104: Coherence selection for States quadrature in F1; implemented by `rho_pos=coherence(spin_system,rho,{{'1H',+1}})`.
- Lines 107-108: tau1 evolution; implemented by `rho_pos=evolution(spin_system,L,[],rho_pos,tau1,1,'final')`.
- Lines 111-112: t1 evolution; implemented by `rho_stack_pos=evolution(spin_system,L,[],rho_pos,t1.timestep/2,t1.nsteps-1,'trajectory')`.
- Lines 115-116: Inversion pulse on 13CA; implemented by `rho_stack_pos=step(spin_system,CAx,rho_stack_pos,pi)`.
- Lines 119-120: t1 rest of the evolution; implemented by `rho_stack_pos=evolution(spin_system,L,[],rho_stack_pos,t1.timestep/2,t1.nsteps-1,'refocus')`.
- Lines 123-124: Inversion pulse on 1H; implemented by `rho_stack_pos=step(spin_system,Hx,rho_stack_pos,pi)`.

### Control flow inferred from the code

- Line 71: conditional branch on `~isfield(parameters,'rho0')`.
- Line 77: conditional branch on `~isfield(parameters,'coil')`.
- Line 215: `for` loop over `name_idx=1:numel(fid_names)`.

### Key state/data transformations

- Lines 54: computes `L` using `L=H+1i*R+1i*K`.
- Lines 57: computes `t1.nsteps` using `t1.nsteps=parameters.npoints(1); t1.timestep=1./parameters.sweep(1)`.
- Lines 58: computes `t2.nsteps` using `t2.nsteps=parameters.npoints(2); t2.timestep=1./parameters.sweep(2)`.
- Lines 59: computes `t3.nsteps` using `t3.nsteps=parameters.npoints(3); t3.timestep=1./parameters.sweep(3)`.
- Lines 62: computes `J_ch` using `J_ch=140; J_nh=92`.
- Lines 63: computes `tau1` using `tau1=abs(1/(4*J_ch))`.
- Lines 64: computes `tau2` using `tau2=abs(1/(4*J_nh))`.
- Lines 65: computes `delta1` using `delta1=abs(1/(4*J_ch))`.
- Lines 66: computes `delta2` using `delta2=12.5e-3`.
- Lines 67: computes `delta3` using `delta3=abs(1/(4*J_nh))`.
- Lines 68: computes `delta4` using `delta4=23.0e-3`.
- Lines 72: computes `HAs` using `HAs=ismember(spin_system.comp.labels,{'HA','HA1','HA2','HA3'})`.
- Lines 73: computes `parameters.rho0` using `parameters.rho0=state(spin_system,'Lz',find(HAs),'cheap')`.
- Lines 78: computes `HNs` using `HNs=ismember(spin_system.comp.labels,{'H'})`.
- Lines 79: computes `parameters.coil` using `parameters.coil=state(spin_system,'L+',find(HNs),'cheap')`.
- Lines 83: computes `Hp` using `Hp=operator(spin_system,'L+',parameters.spins{1})`.
- Lines 84: computes `Hx` using `Hx=(Hp+Hp')/2; Hy=(Hp-Hp')/2i`.
- Lines 87: computes `CAs` using `CAs=strcmp('CA',spin_system.comp.labels)`.

### Local helper functions

- Line 222: `grumble()` — `function grumble(spin_system,parameters,H,R,K)`.
  - Representative operation: `if ~ismember(spin_system.bas.formalism,{'sphten-liouv'})`.
  - Representative operation: `error('this function is only available for sphten-liouv formalism.')`.

## Parameters / inputs

- parameters.npoints -a vector of three integers giving the
- number of points in the three temporal
- dimensions, ordered as [t1 t2 t3].
- parameters.sweep -a vector of three real numbers giving
- the sweep widths in the three frequen-
- cy dimensions, ordered as [f1 f2 f3].
- parameters.spins -isotopes affected by ideal broadband
- pulses, specified as a cell array of
- strings, e.g. {'1H','15N','1H'}.
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

- Protein-specific H(CA)NH experiment (Figure 7.37 of "Protein NMR
- Spectroscopy", 2nd edition) using pre-set values of J-couplings
- used in the magnetisation transfer stages. The simulation uses
- the bidirectional propagation method described in
- The sequence is hard-wired to work on 1H,13C,15N proteins and uses
- PDB labels to select spins that will be affected by otherwise ideal
- pulses. F1 is 1H, F2 is 15N, F3 is 1H. Syntax:
- fid=hcanh(spin_system,parameters,H,R,K)
- parameters.npoints -a vector of three integers giving the
- number of points in the three temporal
- dimensions, ordered as [t1 t2 t3].
- parameters.sweep -a vector of three real numbers giving

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `grumble()`, `isfield()`, `ismember()`, `state()`, `operator()`, `strcmp()`, `step()`, `coherence()`, `evolution()`, `decouple()`, `stitch()`, `fieldnames()`, `ismatrix()`, `all()`, `isvector()`, `any()`.
