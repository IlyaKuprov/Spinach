# experiments/cp_contact_hard.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/experiments/cp_contact_hard.m`
- Signature: `contact_curve=cp_contact_hard(spin_system,parameters,H,R,K)`
- Total lines: 155

## Purpose

Cross-polarisation experiment in the rotating frame. Applies an ideal pi/2 pulse using the specified operators, then evolves the system with the specified spin-lock terms added to the Liovilli- an. The contact curve is returned. Syntax: contact_curve=cp_contact_hard(spin_system,parameters,H,R,K)

## Physical / mathematical content

- This file belongs to the `experiments` part of Spinach. Its role should be read together with nearby files in the same directory, which usually share a common physical regime or infrastructure purpose.

## Numerical / algorithmic content

- Time propagation is explicit. In Spinach this usually means repeated application of matrix exponentials or propagator factorizations to density operators or state vectors in Hilbert/Liouville/Fokker-Planck space.
- The file contains an explicit `grumble(...)` validator, which is Spinach convention for front-loading dimension, type, and regime checks before expensive linear-algebra work begins.
- The file also defines local helper function(s): `grumble()`. This usually means the public entry point is supported by tightly coupled validation or helper logic kept private to the file.

## Code-derived implementation details

### Comment-guided execution stages

- Lines 51-52: Check consistency; implemented by `grumble(parameters,H,R,K)`.
- Lines 54-55: Compose the Liouvillian; implemented by `L=H+1i*R+1i*K`.
- Lines 57-58: Build and project the excitation operator; implemented by `exc_oper=parameters.exc_opers{1}`.
- Lines 64-65: Apply a perfect hard pi/2 excitation pulse; implemented by `rho=step(spin_system,exc_oper,parameters.rho0,pi/2)`.
- Lines 67-69: Preallocate contact curves; implemented by `contact_curve=zeros(size(parameters.coil,2), numel(parameters.time_steps)+1)`.
- Lines 72-73: Loop over the time steps; implemented by `for n=1:numel(parameters.time_steps)`.
- Lines 75-77: Build and project the spin-lock operator; implemented by `irr_oper=2*pi*parameters.irr_powers(1,n)* parameters.irr_opers{1}`.
- Lines 84-85: Apply the evolution step; implemented by `rho=step(spin_system,L+irr_oper,rho,parameters.time_steps(n))`.
- Lines 87-88: Project out the observable; implemented by `contact_curve(:,n+1)=parameters.coil'*rho`.

### Control flow inferred from the code

- Line 59: `for` loop over `n=2:numel(parameters.exc_opers)`.
- Line 73: `for` loop over `n=1:numel(parameters.time_steps)`.
- Line 78: `for` loop over `k=2:numel(parameters.irr_opers)`.

### Key state/data transformations

- Lines 55: computes `L` using `L=H+1i*R+1i*K`.
- Lines 58: computes `exc_oper` using `exc_oper=parameters.exc_opers{1}`.
- Lines 65: computes `rho` using `rho=step(spin_system,exc_oper,parameters.rho0,pi/2)`.
- Lines 68-69: computes `contact_curve` using `contact_curve=zeros(size(parameters.coil,2), numel(parameters.time_steps)+1)`.
- Lines 70: computes `contact_curve(:,1)` using `contact_curve(:,1)=parameters.coil'*rho`.
- Lines 76-77: computes `irr_oper` using `irr_oper=2*pi*parameters.irr_powers(1,n)* parameters.irr_opers{1}`.
- Lines 88: computes `contact_curve(:,n+1)` using `contact_curve(:,n+1)=parameters.coil'*rho`.

### Local helper functions

- Line 95: `grumble()` — `function grumble(parameters,H,R,K)`.
  - Representative operation: `if (~isnumeric(H))||(~isnumeric(R))||(~isnumeric(K))|| (~ismatrix(H))||(~ismatrix(R))||(~ismatrix(K))`.
  - Representative operation: `(~ismatrix(H))||(~ismatrix(R))||(~ismatrix(K))`.

## Parameters / inputs

- parameters.irr_powers -a matrix containing the values
- of the spin-lock nutation fre-
- quency on each channel (rows)
- at each time slice (cols), Hz
- parameters.irr_opers -a cell array of spin operators
- corresponding to the spin-lock
- on each channel
- parameters.exc_opers -a cell array of spin operators
- for the ideal pi/2 excitation
- pulse (same flip angle on all
- channels)
- parameters.time_steps -a vector of time slice durati-
- ons, seconds
- parameters.rho0 -initial state vector
- parameters.coil -detection state vector
- H -Hamiltonian matrix, received from context function
- R -relaxation superoperator, received from context function
- K -kinetics superoperator, received from context function
- Output:
- contact_curve -contact curve detected on the coil
- state specified in parameters.coil

## Implementation structure

- Cross-polarisation experiment in the rotating frame. Applies an
- ideal pi/2 pulse using the specified operators, then evolves the
- system with the specified spin-lock terms added to the Liovilli-
- an. The contact curve is returned. Syntax:
- contact_curve=cp_contact_hard(spin_system,parameters,H,R,K)
- parameters.irr_powers -a matrix containing the values
- of the spin-lock nutation fre-
- quency on each channel (rows)
- at each time slice (cols), Hz
- parameters.irr_opers -a cell array of spin operators
- corresponding to the spin-lock
- on each channel

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `grumble()`, `speye()`, `step()`, `contact_curve()`, `ismatrix()`, `all()`, `isfield()`, `iscell()`, `isscalar()`, `isrow()`, `any()`.
