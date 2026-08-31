# kernel/utilities/expmint2.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/kernel/utilities/expmint2.m`
- Signature: `I=expmint2(spin_system,A,B,C,D,E,T)`
- Total lines: 81

## Purpose

Computes the nested matrix exponential double integral: Integrate[expm(-i*A*(T-t))*B* Integrate[expm(-i*C*(t-x))*D*expm(-i*E*x),{x,0,t}],{t,0,T}] This corresponds to the (1,3) block of the exponential of the auxiliary matrix (http://dx.doi.org/10.1109/TAC.1978.1101743). Syntax: I=expmint2(spin_system,A,B,C,D,E,T)

## Physical / mathematical content

- General mathematical and infrastructure utilities. This area contains finite differences, perturbation theory, graph algorithms, spectral densities, tensor algebra, hash/report helpers, and other reusable numerical components.

## Numerical / algorithmic content

- Time propagation is explicit. In Spinach this usually means repeated application of matrix exponentials or propagator factorizations to density operators or state vectors in Hilbert/Liouville/Fokker-Planck space.
- The file contains an explicit `grumble(...)` validator, which is Spinach convention for front-loading dimension, type, and regime checks before expensive linear-algebra work begins.
- The file also defines local helper function(s): `size()`. This usually means the public entry point is supported by tightly coupled validation or helper logic kept private to the file.

## Code-derived implementation details

### Comment-guided execution stages

- Lines 28-29: Check consistency; implemented by `grumble(A,B,C,D,E,T)`.
- Lines 31-32: Zero filler block; implemented by `Z=sparse(size(A,1),size(A,2))`.
- Lines 34-37: Auxiliary matrix; implemented by `auxmat=[A -1i*B, Z; Z C -1i*D; Z Z E]`.
- Lines 39-40: Exponentiate the auxiliary matrix; implemented by `P=propagator(spin_system,auxmat,T)`.
- Lines 42-43: Build block extractors; implemented by `BE1=[speye(size(A)) Z Z]`.
- Lines 46-47: Extract; implemented by `I=BE1*P*BE3`.

### Key state/data transformations

- Lines 32: computes `Z` using `Z=sparse(size(A,1),size(A,2))`.
- Lines 35-37: computes `auxmat` using `auxmat=[A -1i*B, Z; Z C -1i*D; Z Z E]`.
- Lines 40: computes `P` using `P=propagator(spin_system,auxmat,T)`.
- Lines 43: computes `BE1` using `BE1=[speye(size(A)) Z Z]`.
- Lines 44: computes `BE3` using `BE3=[Z; Z; speye(size(A))]`.
- Lines 47: computes `I` using `I=BE1*P*BE3`.

### Local helper functions

- Line 52: `grumble()` — `function grumble(A,B,C,D,E,T)`.
  - Representative operation: `if (~isnumeric(A))||(~isnumeric(B))||(~isnumeric(C))|| (~isnumeric(D))||(~isnumeric(E))||(~isnumeric(T))`.
  - Representative operation: `(~isnumeric(D))||(~isnumeric(E))||(~isnumeric(T))`.

## Parameters / inputs

- A,B,C,D,E -square matrices
- T -upper limit of the outer integral
- Output:
- I -the integral as above

## Implementation structure

- Computes the nested matrix exponential double integral:
- Integrate[expm(-i*A*(T-t))*B*
- Integrate[expm(-i*C*(t-x))*D*expm(-i*E*x),{x,0,t}],{t,0,T}]
- This corresponds to the (1,3) block of the exponential of the auxiliary
- matrix (http://dx.doi.org/10.1109/TAC.1978.1101743). Syntax:
- I=expmint2(spin_system,A,B,C,D,E,T)
- A,B,C,D,E -square matrices
- T -upper limit of the outer integral
- Output:
- I -the integral as above
- Check consistency
- Zero filler block

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `grumble()`, `propagator()`, `speye()`, `ismatrix()`, `isscalar()`.
