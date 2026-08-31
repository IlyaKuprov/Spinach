# kernel/kinetics/equilibrate.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/kernel/kinetics/equilibrate.m`
- Signature: `c=equilibrate(K,c0)`
- Total lines: 95

## Purpose

Equilibrates linear chemical kinetics and returns a vector of equilibrium concentrations. Syntax: c=equilibrate(K,c0)

## Physical / mathematical content

- This file belongs to the `kernel` part of Spinach. Its role should be read together with nearby files in the same directory, which usually share a common physical regime or infrastructure purpose.

## Numerical / algorithmic content

- The file contains an explicit `grumble(...)` validator, which is Spinach convention for front-loading dimension, type, and regime checks before expensive linear-algebra work begins.
- The file also defines local helper function(s): `size()`. This usually means the public entry point is supported by tightly coupled validation or helper logic kept private to the file.

## Code-derived implementation details

### Comment-guided execution stages

- Lines 24-25: Check consistency; implemented by `grumble(K,c0)`.
- Lines 27-28: Shortcut for zero concentrations; implemented by `if norm(c0,2)==0, c=c0; return; end`.
- Lines 30-31: Recursive calls for independent reactions; implemented by `indep_rxn_idx=scomponents(logical(K)|logical(K)')`.
- Lines 42-43: Assemble the steady state system; implemented by `A=vertcat(ones(1,size(K,2)),K)`.
- Lines 46-47: Check condition number; implemented by `if cond(A)>1/sqrt(eps('double'))`.
- Lines 51-52: Solve the system; implemented by `c=A\b`.

### Control flow inferred from the code

- Line 28: conditional branch on `norm(c0,2)==0, c=c0; return; end`.
- Line 33: conditional branch on `n_indep_rxns>1`.
- Line 35: `for` loop over `n=1:n_indep_rxns`.
- Line 47: conditional branch on `cond(A)>1/sqrt(eps('double'))`.

### Key state/data transformations

- Lines 31: computes `indep_rxn_idx` using `indep_rxn_idx=scomponents(logical(K)|logical(K)')`.
- Lines 32: computes `n_indep_rxns` using `n_indep_rxns=numel(unique(indep_rxn_idx))`.
- Lines 34: computes `c` using `c=zeros(size(c0))`.
- Lines 43: computes `A` using `A=vertcat(ones(1,size(K,2)),K)`.
- Lines 44: computes `b` using `b=vertcat(sum(c0),zeros(size(K,1),1))`.

### Local helper functions

- Line 57: `grumble()` — `function grumble(K,c0)`.
  - Representative operation: `if (~isnumeric(K))||(~isreal(K))||(size(K,1)~=size(K,2))`.
  - Representative operation: `error('K must be a real square matrix.')`.

## Parameters / inputs

- K -reaction rate matrix corresponding to
- dc/dt=K*c, where c is the concentration
- vector
- c0 -vector of initial concentrations

## Outputs

- c -vector of equilibrium concentrations

## Implementation structure

- Equilibrates linear chemical kinetics and returns a vector of
- equilibrium concentrations. Syntax:
- c=equilibrate(K,c0)
- K -reaction rate matrix corresponding to
- dc/dt=K*c, where c is the concentration
- vector
- c0 -vector of initial concentrations
- c -vector of equilibrium concentrations
- Check consistency
- Shortcut for zero concentrations
- Recursive calls for independent reactions
- Assemble the steady state system

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `grumble()`, `scomponents()`, `logical()`, `vertcat()`, `cond()`, `eps()`, `all()`, `iscolumn()`, `any()`.
