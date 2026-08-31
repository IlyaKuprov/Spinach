# kernel/utilities/expmint.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/kernel/utilities/expmint.m`
- Signature: `R=expmint(spin_system,A,B,C,T)`
- Total lines: 91

## Purpose

Computes matrix exponential integrals of the following general type: Integrate[expm(-i*A*t)*B*expm(i*C*t),{t,0,T}] Matrix A must be Hermitian. For further info see the paper by Char- les van Loan (http://dx.doi.org/10.1109/TAC.1978.1101743). Syntax: R=expmint(spin_system,A,B,C,T)

## Physical / mathematical content

- General mathematical and infrastructure utilities. This area contains finite differences, perturbation theory, graph algorithms, spectral densities, tensor algebra, hash/report helpers, and other reusable numerical components.

## Numerical / algorithmic content

- Time propagation is explicit. In Spinach this usually means repeated application of matrix exponentials or propagator factorizations to density operators or state vectors in Hilbert/Liouville/Fokker-Planck space.
- The file contains an explicit `grumble(...)` validator, which is Spinach convention for front-loading dimension, type, and regime checks before expensive linear-algebra work begins.
- The file also defines local helper function(s): `grumble()`. This usually means the public entry point is supported by tightly coupled validation or helper logic kept private to the file.

## Code-derived implementation details

### Comment-guided execution stages

- Lines 34-35: Check consistency; implemented by `grumble(A,B,C,T)`.
- Lines 37-38: Zero argument shortcut; implemented by `if (T==0)||(nnz(B)==0), R=spalloc(size(B,1),size(B,2),0); return; end`.
- Lines 40-41: Build auxiliary matrix; implemented by `auxmat=[-A, 1i*B; 0*A, -C]`.
- Lines 43-44: Get block extraction multipliers; implemented by `BE1=spdiags(ones(size(A,1),1), 0 ,2*size(A,1),size(A,2))`.
- Lines 47-48: Reclaim memory; implemented by `clear('A','B','C')`.
- Lines 50-51: Exponentiate the auxiliary matrix; implemented by `auxmat=propagator(spin_system,auxmat,T)`.
- Lines 53-54: Inform the user; implemented by `report(spin_system,'processing auxiliary matrix blocks ')`.
- Lines 56-57: Cut the blocks out; implemented by `P=(BE1'*auxmat*BE1)'; Q=(BE1'*auxmat*BE2)`.
- Lines 59-60: Reclaim memory; implemented by `clear('auxmat','BE1','BE2')`.
- Lines 62-63: Multiply up the blocks; implemented by `R=clean_up(spin_system,P*Q,spin_system.tols.prop_chop)`.
- Lines 65-66: Reclaim memory; implemented by `clear('P','Q')`.

### Control flow inferred from the code

- Line 38: conditional branch on `(T==0)||(nnz(B)==0), R=spalloc(size(B,1),size(B,2),0); return; end`.

### Key state/data transformations

- Lines 41: computes `auxmat` using `auxmat=[-A, 1i*B; 0*A, -C]`.
- Lines 44: computes `BE1` using `BE1=spdiags(ones(size(A,1),1), 0 ,2*size(A,1),size(A,2))`.
- Lines 45: computes `BE2` using `BE2=spdiags(ones(size(A,1),1),-size(A,1),2*size(A,1),size(A,2))`.
- Lines 57: computes `P` using `P=(BE1'*auxmat*BE1)'; Q=(BE1'*auxmat*BE2)`.
- Lines 63: computes `R` using `R=clean_up(spin_system,P*Q,spin_system.tols.prop_chop)`.

### Local helper functions

- Line 71: `grumble()` — `function grumble(A,B,C,T)`.
  - Representative operation: `if (~isnumeric(A))||(~isnumeric(B))||(~isnumeric(C))|| (~ismatrix(A))||(~ismatrix(B))||(~ismatrix(C))`.
  - Representative operation: `(~ismatrix(A))||(~ismatrix(B))||(~ismatrix(C))`.

## Parameters / inputs

- A,B,C -the three matrices involved in the integral
- T -integration time
- Output:
- R -the resulting integral
- Note: the auxiliary matrix method is massively faster than either
- commutator series or diagonalisation.
- Note: this is the most memory-intensive stage in a lot of calcula-
- tions; memory recycling is aggressive.

## Implementation structure

- Computes matrix exponential integrals of the following general type:
- Integrate[expm(-i*A*t)*B*expm(i*C*t),{t,0,T}]
- Matrix A must be Hermitian. For further info see the paper by Char-
- les van Loan (http://dx.doi.org/10.1109/TAC.1978.1101743). Syntax:
- R=expmint(spin_system,A,B,C,T)
- A,B,C -the three matrices involved in the integral
- T -integration time
- Output:
- R -the resulting integral
- Note: the auxiliary matrix method is massively faster than either
- commutator series or diagonalisation.
- Note: this is the most memory-intensive stage in a lot of calcula-

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `grumble()`, `nnz()`, `spalloc()`, `spdiags()`, `clear()`, `propagator()`, `report()`, `clean_up()`, `ismatrix()`, `all()`, `ishermitian()`, `isscalar()`.
