# experiments/spin_chem/rydmr.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/experiments/spin_chem/rydmr.m`
- Signature: `A=rydmr(spin_system,parameters,H,R,K)`
- Total lines: 77

## Purpose

Singlet-singlet RYDMR experiment using the full kinetics superoper- ator -computes the singlet yield of a radical pair recombination reaction. Syntax: A=rydmr(spin_system,parameters,H,R,K) where H is the Hamiltonian commutation superoperator in zero ex- ternal field, R is the relaxation superoperator and K is the che- mical kinetics superoperator. Parameters: parameters.tol -BICG solver tolerance, 1e-2 is generally g

## Physical / mathematical content

- Spin-chemistry experiment implementations. These routines couple spin evolution to chemical kinetics, radical-pair recombination, exchange, and spin-selective reaction channels.
- The relevant state manifold is the singlet/triplet decomposition, where permutation symmetry controls selection rules, relaxation susceptibility, and convertibility to ordinary magnetisation.

## Numerical / algorithmic content

- The implementation explicitly addresses performance engineering through parallel or GPU execution, which matters because Spinach operators can become extremely large after basis expansion or powder/spatial lifting.
- The file contains an explicit `grumble(...)` validator, which is Spinach convention for front-loading dimension, type, and regime checks before expensive linear-algebra work begins.
- The file also defines local helper function(s): `grumble()`. This usually means the public entry point is supported by tightly coupled validation or helper logic kept private to the file.

## Code-derived implementation details

### Comment-guided execution stages

- Lines 26-27: Check consistency; implemented by `grumble(spin_system,parameters,H,R,K)`.
- Lines 29-31: Get the two-electron singlet state; implemented by `S=singlet(spin_system,spin_system.chem.rp_electrons(1), spin_system.chem.rp_electrons(2))`.
- Lines 33-34: Compose Liouvillian; implemented by `L=H+1i*R+1i*K`.
- Lines 36-37: Normalize the singlet; implemented by `S=S/norm(S,2)`.
- Lines 39-40: Move to GPU if needed; implemented by `if ismember('gpu',spin_system.sys.enable)`.
- Lines 44-46: Compute singlet yield; implemented by `A=spin_system.chem.rp_rates(1)* imag(S'*bicg(L,S,parameters.tol,numel(S)))`.
- Lines 48-49: Gather from GPU if needed; implemented by `if ismember('gpu',spin_system.sys.enable)`.

### Control flow inferred from the code

- Line 40: conditional branch on `ismember('gpu',spin_system.sys.enable)`.
- Line 49: conditional branch on `ismember('gpu',spin_system.sys.enable)`.

### Key state/data transformations

- Lines 30-31: computes `S` using `S=singlet(spin_system,spin_system.chem.rp_electrons(1), spin_system.chem.rp_electrons(2))`.
- Lines 34: computes `L` using `L=H+1i*R+1i*K`.
- Lines 45-46: computes `A` using `A=spin_system.chem.rp_rates(1)* imag(S'*bicg(L,S,parameters.tol,numel(S)))`.

### Local helper functions

- Line 56: `grumble()` — `function grumble(spin_system,parameters,H,R,K)`.
  - Representative operation: `if (~isnumeric(H))||(~isnumeric(R))||(~isnumeric(K))|| (~ismatrix(H))||(~ismatrix(R))||(~ismatrix(K))`.
  - Representative operation: `(~ismatrix(H))||(~ismatrix(R))||(~ismatrix(K))`.

## Outputs

- A -fractional singlet yield

## Implementation structure

- Singlet-singlet RYDMR experiment using the full kinetics superoper-
- ator -computes the singlet yield of a radical pair recombination
- reaction. Syntax:
- A=rydmr(spin_system,parameters,H,R,K)
- where H is the Hamiltonian commutation superoperator in zero ex-
- ternal field, R is the relaxation superoperator and K is the che-
- mical kinetics superoperator. Parameters:
- parameters.tol - BICG solver tolerance,
- 1e-2 is generally good
- A -fractional singlet yield
- Check consistency
- Get the two-electron singlet state

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `grumble()`, `singlet()`, `ismember()`, `gpuArray()`, `bicg()`, `gather()`, `ismatrix()`, `all()`, `isfield()`, `isscalar()`.
