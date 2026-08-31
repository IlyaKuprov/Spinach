# kernel/states/four_spin_states.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/kernel/states/four_spin_states.m`
- Signature: `rho=four_spin_states(spin_system,spins,spin_state)`
- Total lines: 194

## Purpose

Returns user-specified states for a system of four spin-1/2 particles; see also the enclosed Mathematica file. Syntax: rho=four_spin_states(spin_system,spins,spin_state)

## Physical / mathematical content

- State-construction utilities. These routines build equilibrium states, singlets, triplets, partner-state expansions, and physically meaningful density operators in the active basis.
- The relevant state manifold is the singlet/triplet decomposition, where permutation symmetry controls selection rules, relaxation susceptibility, and convertibility to ordinary magnetisation.

## Numerical / algorithmic content

- The file is built around the standard Spinach workflow: create the spin system, choose a basis or context, assemble operators/superoperators, then propagate or analyse the resulting dynamics.
- The file contains an explicit `grumble(...)` validator, which is Spinach convention for front-loading dimension, type, and regime checks before expensive linear-algebra work begins.
- The file also defines local helper function(s): `grumble()`. This usually means the public entry point is supported by tightly coupled validation or helper logic kept private to the file.

## Code-derived implementation details

### Comment-guided execution stages

- Lines 25-26: Check consistency; implemented by `grumble(spin_system,spins,spin_state)`.
- Lines 28-29: Component operators: four-spin; implemented by `XXXX=state(spin_system,{'Lx','Lx','Lx','Lx'},num2cell(spins(:)))`.
- Lines 39-40: Component operators: three-spin; implemented by `XXEZ=state(spin_system,{'Lx','Lx','E' ,'Lz'},num2cell(spins(:)))`.
- Lines 53-54: Component operators: two-spin; implemented by `EEXX=state(spin_system,{'E', 'E' ,'Lx','Lx'},num2cell(spins(:)))`.
- Lines 65-66: Component operators: one-spin; implemented by `EEEZ=state(spin_system,{'E', 'E' ,'E' ,'Lz'},num2cell(spins(:)))`.
- Lines 71-72: Component operators: unit state; implemented by `EEEE=state(spin_system,{'E', 'E', 'E', 'E' },num2cell(spins(:)))`.
- Lines 74-75: Build the state; implemented by `switch spin_state`.
- Lines 163-164: Complain and bomb out; implemented by `error('unknown spin state.')`.

### Control flow inferred from the code

- Line 75: dispatches on `spin_state`; cases `'S(x)S'`, `'S(x)TU'`, `'S(x)T0'`, `'S(x)TD'`, `'TU(x)S'`, `'T0(x)S'`, `'TD(x)S'`, `'TU(x)TU'`, `'T0(x)TU'`, `'TD(x)TU'`, ….

### Key state/data transformations

- Lines 29: computes `XXXX` using `XXXX=state(spin_system,{'Lx','Lx','Lx','Lx'},num2cell(spins(:)))`.
- Lines 30: computes `XXYY` using `XXYY=state(spin_system,{'Lx','Lx','Ly','Ly'},num2cell(spins(:)))`.
- Lines 31: computes `XXZZ` using `XXZZ=state(spin_system,{'Lx','Lx','Lz','Lz'},num2cell(spins(:)))`.
- Lines 32: computes `YYXX` using `YYXX=state(spin_system,{'Ly','Ly','Lx','Lx'},num2cell(spins(:)))`.
- Lines 33: computes `YYYY` using `YYYY=state(spin_system,{'Ly','Ly','Ly','Ly'},num2cell(spins(:)))`.
- Lines 34: computes `YYZZ` using `YYZZ=state(spin_system,{'Ly','Ly','Lz','Lz'},num2cell(spins(:)))`.
- Lines 35: computes `ZZXX` using `ZZXX=state(spin_system,{'Lz','Lz','Lx','Lx'},num2cell(spins(:)))`.
- Lines 36: computes `ZZYY` using `ZZYY=state(spin_system,{'Lz','Lz','Ly','Ly'},num2cell(spins(:)))`.
- Lines 37: computes `ZZZZ` using `ZZZZ=state(spin_system,{'Lz','Lz','Lz','Lz'},num2cell(spins(:)))`.
- Lines 40: computes `XXEZ` using `XXEZ=state(spin_system,{'Lx','Lx','E' ,'Lz'},num2cell(spins(:)))`.
- Lines 41: computes `XXZE` using `XXZE=state(spin_system,{'Lx','Lx','Lz','E' },num2cell(spins(:)))`.
- Lines 42: computes `YYEZ` using `YYEZ=state(spin_system,{'Ly','Ly','E' ,'Lz'},num2cell(spins(:)))`.
- Lines 43: computes `YYZE` using `YYZE=state(spin_system,{'Ly','Ly','Lz','E' },num2cell(spins(:)))`.
- Lines 44: computes `ZZEZ` using `ZZEZ=state(spin_system,{'Lz','Lz','E' ,'Lz'},num2cell(spins(:)))`.
- Lines 45: computes `ZZZE` using `ZZZE=state(spin_system,{'Lz','Lz','Lz','E' },num2cell(spins(:)))`.
- Lines 46: computes `EZXX` using `EZXX=state(spin_system,{'E' ,'Lz','Lx','Lx'},num2cell(spins(:)))`.
- Lines 47: computes `EZYY` using `EZYY=state(spin_system,{'E' ,'Lz','Ly','Ly'},num2cell(spins(:)))`.
- Lines 48: computes `EZZZ` using `EZZZ=state(spin_system,{'E' ,'Lz','Lz','Lz'},num2cell(spins(:)))`.

### Local helper functions

- Line 171: `grumble()` — `function grumble(spin_system,spins,spin_state)`.
  - Representative operation: `if (~isnumeric(spins))||(~isvector(spins))||(numel(spins)~=4)|| (~isreal(spins))||any(spins<1,'all')||any(mod(spins,1)~=0,'all')|| (numel(unique(spins))~=numel(spins))||…`.
  - Representative operation: `(~isreal(spins))||any(spins<1,'all')||any(mod(spins,1)~=0,'all')|| (numel(unique(spins))~=numel(spins))||(~isrow(spins))`.

## Parameters / inputs

- spins -a row vector of four spin numbers
- spin_state -one of the possible singlet-triplet
- product states, see below

## Outputs

- rho -a density matrix (Hilbert space) or
- a state vector (Liouville space)

## Implementation structure

- Returns user-specified states for a system of four spin-1/2
- particles; see also the enclosed Mathematica file. Syntax:
- rho=four_spin_states(spin_system,spins,spin_state)
- spins -a row vector of four spin numbers
- spin_state -one of the possible singlet-triplet
- product states, see below
- rho -a density matrix (Hilbert space) or
- a state vector (Liouville space)
- Check consistency
- Component operators: four-spin
- Component operators: three-spin
- Component operators: two-spin

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `grumble()`, `state()`, `num2cell()`, `spins()`, `isvector()`, `any()`, `isrow()`, `ischar()`.
