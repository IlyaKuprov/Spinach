# experiments/nmr_liquids/noesyhsqc.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/experiments/nmr_liquids/noesyhsqc.m`
- Signature: `fid=noesyhsqc(spin_system,parameters,H,R,K)`
- Total lines: 198

## Purpose

Phase-sensitive NOESY-HSQC sequence described in: The sequence is hard-wired to {F1,F2,F3}={1H,15N,1H} with carbon decoupled throughout. Syntax: fid=noesyhsqc(spin_system,parameters,H,R,K)

## Physical / mathematical content

- Liquid-state pulse sequence implementations. These are production experiment kernels that carry out coherence transfer, mixing, refocusing, decoupling, and indirect evolution on precomputed Hamiltonian/relaxation/kinetics operators.
- Orientation or trajectory averaging is performed numerically, so grid design, weights, and integration error control matter directly to accuracy and runtime.

## Numerical / algorithmic content

- Time propagation is explicit. In Spinach this usually means repeated application of matrix exponentials or propagator factorizations to density operators or state vectors in Hilbert/Liouville/Fokker-Planck space.
- Numerical integration over angles or geometry is part of the implementation, so point placement and weights are as important as the local Hamiltonian calculations.
- The file contains an explicit `grumble(...)` validator, which is Spinach convention for front-loading dimension, type, and regime checks before expensive linear-algebra work begins.
- The file also defines local helper function(s): `grumble()`. This usually means the public entry point is supported by tightly coupled validation or helper logic kept private to the file.

## Code-derived implementation details

### Comment-guided execution stages

- Lines 53-54: Check consistency; implemented by `grumble(spin_system,parameters,H,R,K)`.
- Lines 56-57: Compose Liouvillian; implemented by `L=H+1i*R+1i*K`.
- Lines 59-60: Decouple carbon; implemented by `if ismember('13C',spin_system.comp.isotopes)`.
- Lines 66-67: Evolution time discretization; implemented by `t1.nsteps=parameters.npoints(1); t1.timestep=1/parameters.sweep(1)`.
- Lines 71-72: J-coupling transfer time; implemented by `delta=abs(1/(4*parameters.J))`.
- Lines 74-75: Initial and detection states; implemented by `rho=state(spin_system,'Lz','1H','cheap')`.
- Lines 78-79: Pulse operators; implemented by `Hx=operator(spin_system,'Lx','1H')`.
- Lines 83-84: Decouple nitrogen during NOESY; implemented by `L_dec=decouple(spin_system,L,[],{'15N'})`.
- Lines 88-89: Run NOESY forward; implemented by `rho=step(spin_system,Hx,rho,pi/2)`.
- Lines 104-105: Run HSQC preparation period forward; implemented by `rho_stack_pos=evolution(spin_system,L,[],rho_stack_pos,delta,1,'final')`.
- Lines 116-117: Run the rest of HSQC backward; implemented by `L_dec=decouple(spin_system,L,[],{'15N'})`.
- Lines 127-128: Stitch the halves; implemented by `report(spin_system,'stitching forward and backward trajectories ')`.
- Lines 134-135: Permute dimensions; implemented by `fid_names=fieldnames(fid)`.

### Control flow inferred from the code

- Line 60: conditional branch on `ismember('13C',spin_system.comp.isotopes)`.
- Line 136: `for` loop over `name_idx=1:numel(fid_names)`.

### Key state/data transformations

- Lines 57: computes `L` using `L=H+1i*R+1i*K`.
- Lines 62: computes `R` using `R=decouple(spin_system,R,[],{'13C'})`.
- Lines 63: computes `K` using `K=decouple(spin_system,K,[],{'13C'})`.
- Lines 67: computes `t1.nsteps` using `t1.nsteps=parameters.npoints(1); t1.timestep=1/parameters.sweep(1)`.
- Lines 68: computes `t2.nsteps` using `t2.nsteps=parameters.npoints(2); t2.timestep=1/parameters.sweep(2)`.
- Lines 69: computes `t3.nsteps` using `t3.nsteps=parameters.npoints(3); t3.timestep=1/parameters.sweep(3)`.
- Lines 72: computes `delta` using `delta=abs(1/(4*parameters.J))`.
- Lines 75: computes `rho` using `rho=state(spin_system,'Lz','1H','cheap')`.
- Lines 76: computes `coil` using `coil=state(spin_system,'L+','1H','cheap')`.
- Lines 79: computes `Hx` using `Hx=operator(spin_system,'Lx','1H')`.
- Lines 80: computes `Hy` using `Hy=operator(spin_system,'Ly','1H')`.
- Lines 81: computes `Nx` using `Nx=operator(spin_system,'Lx','15N')`.
- Lines 84: computes `L_dec` using `L_dec=decouple(spin_system,L,[],{'15N'})`.
- Lines 85: computes `R_dec` using `R_dec=decouple(spin_system,R,[],{'15N'})`.
- Lines 86: computes `K_dec` using `K_dec=decouple(spin_system,K,[],{'15N'})`.
- Lines 90: computes `rho_stack` using `rho_stack=evolution(spin_system,L_dec,[],rho,t1.timestep,t1.nsteps-1,'trajectory')`.
- Lines 91: computes `rho_stack_pos` using `rho_stack_pos=coherence(spin_system,rho_stack,{{'1H',+1}})`.
- Lines 92: computes `rho_stack_neg` using `rho_stack_neg=coherence(spin_system,rho_stack,{{'1H',-1}})`.

### Local helper functions

- Line 143: `grumble()` — `function grumble(spin_system,parameters,H,R,K)`.
  - Representative operation: `if ~ismember(spin_system.bas.formalism,{'sphten-liouv'})`.
  - Representative operation: `error('this function is only available for sphten-liouv formalism.')`.

## Parameters / inputs

- parameters.npoints -a vector of three integers giving the
- number of points in the three temporal
- dimensions, ordered as [t1 t2 t3].
- parameters.sweep -a vector of three real numbers giving
- the sweep widths in the three frequen-
- cy dimensions, ordered as [f1 f2 f3].
- parameters.J -J-coupling for the HSQC stage magneti-
- sation transfer, Hz.
- parameters.tmix -NOESY stage mixing time, seconds.
- H -Hamiltonian matrix, received from context function
- R -relaxation superoperator, received from context function
- K -kinetics superoperator, received from context function

## Outputs

- fid -a structure with four fields: fid.pos_pos, fid.pos_neg,
- fid.neg_pos, fid.neg_neg that are used in the subsequ-
- ent States quadrature processing
- Notes: the sequence starts with a pure Lz on protons at the mo-
- ment and assumes that the relaxation superoperator is not
- thermalised -the relaxation destination is the zero state.
- The first citation is the NOESY-HMQC precursor to this
- NOESY-HSQC implementation.

## Implementation structure

- Phase-sensitive NOESY-HSQC sequence described in:
- The sequence is hard-wired to {F1,F2,F3}={1H,15N,1H} with carbon
- decoupled throughout. Syntax:
- fid=noesyhsqc(spin_system,parameters,H,R,K)
- parameters.npoints -a vector of three integers giving the
- number of points in the three temporal
- dimensions, ordered as [t1 t2 t3].
- parameters.sweep -a vector of three real numbers giving
- the sweep widths in the three frequen-
- cy dimensions, ordered as [f1 f2 f3].
- parameters.J -J-coupling for the HSQC stage magneti-
- sation transfer, Hz.

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `grumble()`, `ismember()`, `decouple()`, `state()`, `operator()`, `step()`, `evolution()`, `coherence()`, `homospoil()`, `report()`, `stitch()`, `fieldnames()`, `ismatrix()`, `all()`, `isfield()`, `isvector()`.
