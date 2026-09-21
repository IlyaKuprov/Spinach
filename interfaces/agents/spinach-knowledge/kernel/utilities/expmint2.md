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
