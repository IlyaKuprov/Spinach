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
