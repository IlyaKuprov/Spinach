# experiments/esr_dipolar/ridme.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/experiments/esr_dipolar/ridme.m`
- Signature: `answer=ridme(spin_system,parameters,H,R,K)`
- Total lines: 182

## Purpose

RIDME pulse sequence. Idealized hard pulses are used, the pulses only affect the user-specified electron. Syntax: answer=ridme(spin_system,parameters,H,R,K) where H is the Hamiltonian commutation superoperator, R is the relaxa- tion superoperator and K is the chemical kinetics superoperator.

## Physical / mathematical content

- Dipolar ESR experiment implementations. The pulse logic resolves dipolar couplings by echo modulation, with selective excitation and time-domain accumulation.
- Orientation or trajectory averaging is performed numerically, so grid design, weights, and integration error control matter directly to accuracy and runtime.
- The spin physics includes through-space magnetic dipole-dipole coupling, a rank-2 anisotropic interaction with strong orientation dependence and characteristic secular/non-secular structure.

## Numerical / algorithmic content

- Time propagation is explicit. In Spinach this usually means repeated application of matrix exponentials or propagator factorizations to density operators or state vectors in Hilbert/Liouville/Fokker-Planck space.
- Numerical integration over angles or geometry is part of the implementation, so point placement and weights are as important as the local Hamiltonian calculations.
- The file contains an explicit `grumble(...)` validator, which is Spinach convention for front-loading dimension, type, and regime checks before expensive linear-algebra work begins.
- The file also defines local helper function(s): `grumble()`. This usually means the public entry point is supported by tightly coupled validation or helper logic kept private to the file.

## Code-derived implementation details

### Comment-guided execution stages

- Lines 53-54: Check consistency; implemented by `grumble(parameters,H,R,K)`.
- Lines 56-57: Compose Liouvillian; implemented by `L=H+1i*R+1i*K`.
- Lines 59-61: Get operators and states; implemented by `coil_imag=(state(spin_system,{'L+'},{parameters.probe_spin})+ state(spin_system,{'L-'},{parameters.probe_spin}))/2`.
- Lines 69-70: First pulse; implemented by `rho=step(spin_system,Sx,parameters.rho0,+pi/2)`.
- Lines 72-73: Evolution (tau 1); implemented by `rho_stack=evolution(spin_system,L,[],rho,(parameters.stepsize*parameters.nsteps(1)),1,'final')`.
- Lines 75-76: Second pulse; implemented by `rho_stack=step(spin_system,Sx,rho_stack,+pi)`.
- Lines 78-80: Evolution (initial tau 1 and tau 2); implemented by `rho_stack=evolution(spin_system,L,[],rho_stack,parameters.stepsize, (parameters.nsteps(1)+parameters.nsteps(2)),'trajectory')`.
- Lines 82-83: Third pulse and a phase cycle; implemented by `rho_stack_px=step(spin_system,Sx,rho_stack,+pi/2)`.
- Lines 88-89: Evolution (mixing time); implemented by `rho_stack_px=evolution(spin_system,L,[],rho_stack_px,parameters.tmix,1,'final')`.
- Lines 94-95: Fourth pulse and a phase cycle; implemented by `rho_stack_pxpx=step(spin_system,Sx,rho_stack_px,+pi/2)`.
- Lines 100-102: Evolution (remaining tau 1 + tau 2); implemented by `rho_stack_pxpx(:,end:-1:1)=evolution(spin_system,L,[],rho_stack_pxpx(:,end:-1:1),parameters.stepsize, (parameters.nsteps(1)+parameters.nsteps(2)),'refocus')`.
- Lines 110-111: Fifth pulse; implemented by `rho_stack_pxpxpx=step(spin_system,Sx,rho_stack_pxpx,+pi)`.
- Lines 116-117: Evolution (tau 2); implemented by `rho_stack_pxpxpx=evolution(spin_system,L,[],rho_stack_pxpxpx,(parameters.stepsize*parameters.nsteps(2)),1,'final')`.
- Lines 122-123: Observation; implemented by `answer.pxpxpx.real=coil_real'*rho_stack_pxpxpx/norm(coil_real,2)`.

### Key state/data transformations

- Lines 57: computes `L` using `L=H+1i*R+1i*K`.
- Lines 60-61: computes `coil_imag` using `coil_imag=(state(spin_system,{'L+'},{parameters.probe_spin})+ state(spin_system,{'L-'},{parameters.probe_spin}))/2`.
- Lines 62-63: computes `coil_real` using `coil_real=(state(spin_system,{'L+'},{parameters.probe_spin})- state(spin_system,{'L-'},{parameters.probe_spin}))/2i`.
- Lines 64-65: computes `Sx` using `Sx=(operator(spin_system,{'L+'},{parameters.probe_spin})+ operator(spin_system,{'L-'},{parameters.probe_spin}))/2`.
- Lines 66-67: computes `Sy` using `Sy=(operator(spin_system,{'L+'},{parameters.probe_spin})- operator(spin_system,{'L-'},{parameters.probe_spin}))/2i`.
- Lines 70: computes `rho` using `rho=step(spin_system,Sx,parameters.rho0,+pi/2)`.
- Lines 73: computes `rho_stack` using `rho_stack=evolution(spin_system,L,[],rho,(parameters.stepsize*parameters.nsteps(1)),1,'final')`.
- Lines 83: computes `rho_stack_px` using `rho_stack_px=step(spin_system,Sx,rho_stack,+pi/2)`.
- Lines 84: computes `rho_stack_py` using `rho_stack_py=step(spin_system,Sy,rho_stack,+pi/2)`.
- Lines 85: computes `rho_stack_mx` using `rho_stack_mx=step(spin_system,Sx,rho_stack,-pi/2)`.
- Lines 86: computes `rho_stack_my` using `rho_stack_my=step(spin_system,Sy,rho_stack,-pi/2)`.
- Lines 95: computes `rho_stack_pxpx` using `rho_stack_pxpx=step(spin_system,Sx,rho_stack_px,+pi/2)`.
- Lines 96: computes `rho_stack_pypy` using `rho_stack_pypy=step(spin_system,Sy,rho_stack_py,+pi/2)`.
- Lines 97: computes `rho_stack_mxmx` using `rho_stack_mxmx=step(spin_system,Sx,rho_stack_mx,-pi/2)`.
- Lines 98: computes `rho_stack_mymy` using `rho_stack_mymy=step(spin_system,Sy,rho_stack_my,-pi/2)`.
- Lines 101-102: computes `rho_stack_pxpx(:,end:-1:1)` using `rho_stack_pxpx(:,end:-1:1)=evolution(spin_system,L,[],rho_stack_pxpx(:,end:-1:1),parameters.stepsize, (parameters.nsteps(1)+parameters.nsteps(2)),'refocus')`.
- Lines 103-104: computes `rho_stack_pypy(:,end:-1:1)` using `rho_stack_pypy(:,end:-1:1)=evolution(spin_system,L,[],rho_stack_pypy(:,end:-1:1),parameters.stepsize, (parameters.nsteps(1)+parameters.nsteps(2)),'refocus')`.
- Lines 105-106: computes `rho_stack_mxmx(:,end:-1:1)` using `rho_stack_mxmx(:,end:-1:1)=evolution(spin_system,L,[],rho_stack_mxmx(:,end:-1:1),parameters.stepsize, (parameters.nsteps(1)+parameters.nsteps(2)),'refocus')`.

### Local helper functions

- Line 135: `grumble()` — `function grumble(parameters,H,R,K)`.
  - Representative operation: `if (~isnumeric(H))||(~isnumeric(R))||(~isnumeric(K))|| (~ismatrix(H))||(~ismatrix(R))||(~ismatrix(K))`.
  - Representative operation: `(~ismatrix(H))||(~ismatrix(R))||(~ismatrix(K))`.

## Parameters / inputs

- H Hamiltonian (received from the context
- function)
- R relaxation superoperator (received from
- the context function)
- K kinetics superoperator (received from
- the context function)
- parameters.rho0 initial state
- parameters.probe_spin number of the spin on which the
- sequence operates
- parameters.stepsize step size for the increment of
- the relaxation period, seconds
- parameters.nsteps(1) number of steps for tau 1
- parameters.nsteps(2) number of steps for tau 2
- parameters.tmix mixing time, seconds

## Outputs

- answer.pxpxpx.(real,imag)
- answer.pypypx.(real,imag)
- answer.mxmxpx.(real,imag)
- answer.mymypx.(real,imag) -quadrature components of the signal
- corresponding to the phase cycle in-
- stances on third, fourth, and fifth
- pulse in the RIDME sequence
- Notes: for this experiment to work, relaxation must be present.

## Implementation structure

- RIDME pulse sequence. Idealized hard pulses are used, the pulses only
- affect the user-specified electron. Syntax:
- answer=ridme(spin_system,parameters,H,R,K)
- where H is the Hamiltonian commutation superoperator, R is the relaxa-
- tion superoperator and K is the chemical kinetics superoperator.
- H Hamiltonian (received from the context
- function)
- R relaxation superoperator (received from
- the context function)
- K kinetics superoperator (received from
- parameters.rho0 initial state
- parameters.probe_spin number of the spin on which the

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `grumble()`, `state()`, `operator()`, `step()`, `evolution()`, `rho_stack_pxpx()`, `rho_stack_pypy()`, `rho_stack_mxmx()`, `rho_stack_mymy()`, `ismatrix()`, `all()`, `isfield()`, `isscalar()`, `isrow()`, `any()`.
