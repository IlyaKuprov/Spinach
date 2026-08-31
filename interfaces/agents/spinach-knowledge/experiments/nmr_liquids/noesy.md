# experiments/nmr_liquids/noesy.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/experiments/nmr_liquids/noesy.m`
- Signature: `fid=noesy(spin_system,parameters,H,R,K)`
- Total lines: 206

## Purpose

Phase-sensitive homonuclear NOESY pulse sequence from:

## Physical / mathematical content

- Liquid-state pulse sequence implementations. These are production experiment kernels that carry out coherence transfer, mixing, refocusing, decoupling, and indirect evolution on precomputed Hamiltonian/relaxation/kinetics operators.

## Numerical / algorithmic content

- Time propagation is explicit. In Spinach this usually means repeated application of matrix exponentials or propagator factorizations to density operators or state vectors in Hilbert/Liouville/Fokker-Planck space.
- The file contains an explicit `grumble(...)` validator, which is Spinach convention for front-loading dimension, type, and regime checks before expensive linear-algebra work begins.
- The file also defines local helper function(s): `grumble()`. This usually means the public entry point is supported by tightly coupled validation or helper logic kept private to the file.

## Code-derived implementation details

### Comment-guided execution stages

- Lines 59-60: Consistency check; implemented by `grumble(spin_system,parameters,H,R,K)`.
- Lines 62-63: Default version includes homospoil; implemented by `if ~isfield(parameters,'oldschool')`.
- Lines 67-68: Coherent evolution timestep; implemented by `timestep=1./parameters.sweep`.
- Lines 70-71: Detection state; implemented by `coil=state(spin_system,'L+',parameters.spins{1},'cheap')`.
- Lines 73-74: Pulse operators; implemented by `Lx=operator(spin_system,'Lx',parameters.spins{1})`.
- Lines 77-78: Decoupling; implemented by `if isfield(parameters,'decouple')`.
- Lines 83-84: First pulse; implemented by `rho=step(spin_system,Lx,parameters.rho0,pi/2)`.
- Lines 86-87: Phase cycle specification; implemented by `Op2={Lx,Ly,Lx,Ly}; An2={+pi/2,+pi/2,-pi/2,-pi/2}`.
- Lines 90-91: FID phase cycle; implemented by `fids=cell(1,4)`.
- Lines 93-94: Phase cycle loop; implemented by `for n=1:4`.
- Lines 96-98: F1 evolution; implemented by `rho_stack=evolution(spin_system,H+1i*R+1i*K,[],rho,timestep(1), parameters.npoints(1)-1,'trajectory')`.
- Lines 99-100: Second pulse; implemented by `rho_stack=step(spin_system,Op2{n},rho_stack,An2{n})`.
- Lines 102-103: Gradient purification; implemented by `if parameters.oldschool`.
- Lines 105-108: Mixing time under full evolution generator, thus including ZQ artefacts and other mess; implemented by `rho_stack=evolution(spin_system,H+1i*R+1i*K,[], rho_stack,parameters.tmix,1,'final')`.
- Lines 112-113: Destroy everything except longitudinal magnetisation; implemented by `rho_stack=homospoil(spin_system,rho_stack,'destroy')`.
- Lines 115-117: Mixing time under relaxation and kinetics; implemented by `rho_stack=evolution(spin_system,1i*R+1i*K,[], rho_stack,parameters.tmix,1,'final')`.
- Lines 121-122: Third pulse; implemented by `rho_stack=step(spin_system,Op3{n},rho_stack,An3{n})`.
- Lines 124-126: F2 evolution and detection; implemented by `fids{n}=evolution(spin_system,H+1i*R+1i*K,coil,rho_stack, timestep(2),parameters.npoints(2)-1,'observable')`.

### Control flow inferred from the code

- Line 63: conditional branch on `~isfield(parameters,'oldschool')`.
- Line 78: conditional branch on `isfield(parameters,'decouple')`.
- Line 94: `for` loop over `n=1:4`.
- Line 103: conditional branch on `parameters.oldschool`.

### Key state/data transformations

- Lines 64: computes `parameters.oldschool` using `parameters.oldschool=false()`.
- Lines 68: computes `timestep` using `timestep=1./parameters.sweep`.
- Lines 71: computes `coil` using `coil=state(spin_system,'L+',parameters.spins{1},'cheap')`.
- Lines 74: computes `Lx` using `Lx=operator(spin_system,'Lx',parameters.spins{1})`.
- Lines 75: computes `Ly` using `Ly=operator(spin_system,'Ly',parameters.spins{1})`.
- Lines 79-80: computes `[H,parameters.rho0]` using `[H,parameters.rho0]=decouple(spin_system,H,parameters.rho0, parameters.decouple)`.
- Lines 84: computes `rho` using `rho=step(spin_system,Lx,parameters.rho0,pi/2)`.
- Lines 87: computes `Op2` using `Op2={Lx,Ly,Lx,Ly}; An2={+pi/2,+pi/2,-pi/2,-pi/2}`.
- Lines 88: computes `Op3` using `Op3={Ly,Ly,Ly,Ly}; An3={+pi/2,+pi/2,+pi/2,+pi/2}`.
- Lines 91: computes `fids` using `fids=cell(1,4)`.
- Lines 97-98: computes `rho_stack` using `rho_stack=evolution(spin_system,H+1i*R+1i*K,[],rho,timestep(1), parameters.npoints(1)-1,'trajectory')`.
- Lines 125-126: computes `fids{n}` using `fids{n}=evolution(spin_system,H+1i*R+1i*K,coil,rho_stack, timestep(2),parameters.npoints(2)-1,'observable')`.
- Lines 131: computes `fid.cos` using `fid.cos=fids{1}-fids{3}; fid.sin=fids{2}-fids{4}`.

### Local helper functions

- Line 136: `grumble()` — `function grumble(spin_system,parameters,H,R,K)`.
  - Representative operation: `if ~ismember(spin_system.bas.formalism,{'sphten-liouv','zeeman-liouv'})`.
  - Representative operation: `error('this function is only available for sphten-liouv and zeeman-liouv formalisms.')`.

## Syntax

```matlab
fid=noesy(spin_system,parameters,H,R,K)
```

## Parameters / inputs

- parameters.sweep -sweep widths, Hz
- parameters.npoints -number of points for both dimensions
- parameters.spins -nuclei on which the sequence runs,
- specified as {'1H'}, {'13C'}, etc.
- parameters.tmix -mixing time, seconds
- parameters.decouple -spins to be decoupled, specified either
- by name, e.g. {'13C','1H'}, or by a list
- of numbers, e.g. [1 2]
- parameters.rho0 -initial state; skip this and specify
- parameters.needs={'rho_eq'} to start
- from exact thermal equilibrium
- parameters.oldschool -set to 1 to disable homospoil gradient
- before the mixing time
- H -Hamiltonian matrix, received from context function
- R -relaxation superoperator, received from context function
- K -kinetics superoperator, received from context function

## Outputs

- fid.cos,fid.sin -two components of the FID for F1 hyper-
- complex processing
- Note: this function is used for extreme simulations (proteins
- and nucleic acids) -its layout is optimised for minimum
- memory footprint rather than CPU time.
- Note: non-empty analytical decoupling is meaningful only in
- sphten-liouv formalism.

## Implementation structure

- Phase-sensitive homonuclear NOESY pulse sequence from:
- fid=noesy(spin_system,parameters,H,R,K)
- parameters.sweep -sweep widths, Hz
- parameters.npoints -number of points for both dimensions
- parameters.spins -nuclei on which the sequence runs,
- specified as {'1H'}, {'13C'}, etc.
- parameters.tmix -mixing time, seconds
- parameters.decouple -spins to be decoupled, specified either
- by name, e.g. {'13C','1H'}, or by a list
- of numbers, e.g. [1 2]
- parameters.rho0 -initial state; skip this and specify
- parameters.needs={'rho_eq'} to start

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `grumble()`, `isfield()`, `false()`, `state()`, `operator()`, `decouple()`, `step()`, `evolution()`, `timestep()`, `homospoil()`, `ismember()`, `ismatrix()`, `all()`, `elseif()`, `any()`, `iscell()`.
