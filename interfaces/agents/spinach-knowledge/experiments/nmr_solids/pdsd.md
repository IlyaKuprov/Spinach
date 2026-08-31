# experiments/nmr_solids/pdsd.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/experiments/nmr_solids/pdsd.m`
- Signature: `fid=pdsd(spin_system,parameters,H,R,K)`
- Total lines: 145

## Purpose

A simplified model of the PDSD experiment using NOESY type quadrature detection and phase cycle. To be cal- led from the singlerot context. Syntax: fid=pdsd(spin_system,parameters,H,R,K)

## Physical / mathematical content

- Solid-state pulse sequence implementations. The core ingredients are anisotropic Hamiltonians, rotor synchronisation, cross-polarisation, recoupling/decoupling, and powder or rotor-stack propagation.
- Orientation or trajectory averaging is performed numerically, so grid design, weights, and integration error control matter directly to accuracy and runtime.

## Numerical / algorithmic content

- Time propagation is explicit. In Spinach this usually means repeated application of matrix exponentials or propagator factorizations to density operators or state vectors in Hilbert/Liouville/Fokker-Planck space.
- Numerical integration over angles or geometry is part of the implementation, so point placement and weights are as important as the local Hamiltonian calculations.
- The file contains an explicit `grumble(...)` validator, which is Spinach convention for front-loading dimension, type, and regime checks before expensive linear-algebra work begins.
- The file also defines local helper function(s): `size()`. This usually means the public entry point is supported by tightly coupled validation or helper logic kept private to the file.

## Code-derived implementation details

### Comment-guided execution stages

- Lines 42-43: Check consistency; implemented by `grumble(parameters,H,R,K)`.
- Lines 45-46: Compose the Liouvillian and its decoupled version; implemented by `L=H+1i*R+1i*K; LD=decouple(spin_system,L,[],{'1H'})`.
- Lines 48-49: Get evolution timesteps; implemented by `timestep=1/parameters.sweep`.
- Lines 51-52: Get pulse operators; implemented by `Cx=operator(spin_system,'Lx','13C')`.
- Lines 59-60: Quadrature detection state; implemented by `coil=state(spin_system,'L+','13C','cheap')`.
- Lines 64-65: Skip CP and start with Ly on carbons; implemented by `rho=state(spin_system,'Ly','13C','cheap')`.
- Lines 68-69: Phase cycle specification; implemented by `Op2={Cx,Cy,Cx,Cy}; An2={+pi/2,+pi/2,-pi/2,-pi/2}`.
- Lines 72-73: FID phase cycle; implemented by `fids=cell(1,4)`.
- Lines 75-76: Phase cycle loop; implemented by `for n=1:4`.
- Lines 78-80: F1 evolution under proton decoupling; implemented by `rho_stack=evolution(spin_system,LD,[],rho,timestep, parameters.npoints(1)-1,'trajectory')`.
- Lines 81-82: Second pulse; implemented by `rho_stack=step(spin_system,Op2{n},rho_stack,An2{n})`.
- Lines 84-86: Mixing time evolution with proton irradiation; implemented by `rho_stack=evolution(spin_system,L+2*pi*parameters.rate*Hx, [],rho_stack,parameters.tmix,1,'final')`.
- Lines 88-89: Third pulse; implemented by `rho_stack=step(spin_system,Op3{n},rho_stack,An3{n})`.
- Lines 91-92: Wipe the proton subspace (decoupling); implemented by `[~,rho_stack]=decouple(spin_system,[],rho_stack,{'1H'})`.
- Lines 94-96: F2 evolution and detection under proton decoupling; implemented by `fids{n}=evolution(spin_system,LD,coil,rho_stack,timestep, parameters.npoints(2)-1,'observable')`.
- Lines 100-101: Axial peak elimination; implemented by `fid.cos=fids{1}-fids{3}; fid.sin=fids{2}-fids{4}`.

### Control flow inferred from the code

- Line 76: `for` loop over `n=1:4`.

### Key state/data transformations

- Lines 46: computes `L` using `L=H+1i*R+1i*K; LD=decouple(spin_system,L,[],{'1H'})`.
- Lines 49: computes `timestep` using `timestep=1/parameters.sweep`.
- Lines 52: computes `Cx` using `Cx=operator(spin_system,'Lx','13C')`.
- Lines 53: computes `Cy` using `Cy=operator(spin_system,'Ly','13C')`.
- Lines 56: computes `Hx` using `Hx=operator(spin_system,'Lx','1H')`.
- Lines 60: computes `coil` using `coil=state(spin_system,'L+','13C','cheap')`.
- Lines 65: computes `rho` using `rho=state(spin_system,'Ly','13C','cheap')`.
- Lines 69: computes `Op2` using `Op2={Cx,Cy,Cx,Cy}; An2={+pi/2,+pi/2,-pi/2,-pi/2}`.
- Lines 70: computes `Op3` using `Op3={Cx,Cx,Cx,Cx}; An3={+pi/2,+pi/2,+pi/2,+pi/2}`.
- Lines 73: computes `fids` using `fids=cell(1,4)`.
- Lines 79-80: computes `rho_stack` using `rho_stack=evolution(spin_system,LD,[],rho,timestep, parameters.npoints(1)-1,'trajectory')`.
- Lines 92: computes `[~,rho_stack]` using `[~,rho_stack]=decouple(spin_system,[],rho_stack,{'1H'})`.
- Lines 95-96: computes `fids{n}` using `fids{n}=evolution(spin_system,LD,coil,rho_stack,timestep, parameters.npoints(2)-1,'observable')`.
- Lines 101: computes `fid.cos` using `fid.cos=fids{1}-fids{3}; fid.sin=fids{2}-fids{4}`.

### Local helper functions

- Line 106: `grumble()` — `function grumble(parameters,H,R,K)`.
  - Representative operation: `if (~isnumeric(H))||(~isnumeric(R))||(~isnumeric(K))|| (~isequal(size(H),size(R),size(K)))|| (size(H,1)~=size(H,2))`.
  - Representative operation: `(~isequal(size(H),size(R),size(K)))|| (size(H,1)~=size(H,2))`.

## Parameters / inputs

- spin_system -Spinach spin system object
- parameters.sweep -sweep width in Hz
- parameters.npoints -two-element vector giving the
- number of complex points in the
- indirect and direct dimensions
- parameters.tmix -mixing time in seconds
- parameters.rate -MAS rate in Hz, used to set
- proton irradiation power
- parameters.spc_dim -spatial dimension of the MAS
- problem, received from the
- context function
- H, R, K -Hamiltonian, relaxation, and
- kinetics superoperators, recei-
- ved from the context function

## Outputs

- fid.cos, fid.sin -States quadrature components
- of the 2D PDSD spectrum

## Implementation structure

- A simplified model of the PDSD experiment using NOESY
- type quadrature detection and phase cycle. To be cal-
- led from the singlerot context. Syntax:
- fid=pdsd(spin_system,parameters,H,R,K)
- spin_system -Spinach spin system object
- parameters.sweep -sweep width in Hz
- parameters.npoints -two-element vector giving the
- number of complex points in the
- indirect and direct dimensions
- parameters.tmix -mixing time in seconds
- parameters.rate -MAS rate in Hz, used to set
- proton irradiation power

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `grumble()`, `decouple()`, `operator()`, `speye()`, `state()`, `evolution()`, `step()`, `isequal()`, `isfield()`, `isscalar()`, `any()`.
