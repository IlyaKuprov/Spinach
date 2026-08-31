# kernel/operators/irr_sph_ten.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/kernel/operators/irr_sph_ten.m`
- Signature: `T=irr_sph_ten(mult,k)`
- Total lines: 108

## Purpose

Single-spin irreducible spherical tensor operators T(k,m) obeying the following commutation relation: [Lz,T_km]=m*T_km

## Physical / mathematical content

- Operator-construction utilities. They build bases and irreducible tensor representations for spin, bosonic, and transition operators.

## Numerical / algorithmic content

- The file contains an explicit `grumble(...)` validator, which is Spinach convention for front-loading dimension, type, and regime checks before expensive linear-algebra work begins.
- The file also defines local helper function(s): `grumble()`. This usually means the public entry point is supported by tightly coupled validation or helper logic kept private to the file.

## Code-derived implementation details

### Comment-guided execution stages

- Lines 37-38: Adapt to the input style; implemented by `if nargin==1`.
- Lines 40-41: Check consistency; implemented by `grumble(mult)`.
- Lines 43-44: Generate tensors of all ranks and put them into a cell array; implemented by `T=cell(mult^2,1)`.
- Lines 46-47: Fill in recursively; implemented by `for n=0:(mult-1)`.
- Lines 53-54: Check consistency; implemented by `grumble(mult,k)`.
- Lines 56-57: Process zero-rank requests; implemented by `if k==0`.
- Lines 59-60: Return a unit matrix; implemented by `T={speye(mult)}`.
- Lines 64-65: Get Pauli matrices; implemented by `L=pauli(mult)`.
- Lines 67-68: Preallocate the cell array; implemented by `T=cell(2*k+1,1)`.
- Lines 70-71: Get the top state; implemented by `T{1}=((-1)^k)*(2^(-k/2))*L.p^k`.
- Lines 73-74: Apply sequential lowering using Racah's commutation rule; implemented by `for n=2:(2*k+1)`.
- Lines 82-83: Catch incorrect calls; implemented by `error('the number of input arguments must be one or two.')`.

### Control flow inferred from the code

- Line 38: conditional branch on `nargin==1`.
- Line 47: `for` loop over `n=0:(mult-1)`.
- Line 57: conditional branch on `k==0`.
- Line 74: `for` loop over `n=2:(2*k+1)`.

### Key state/data transformations

- Lines 44: computes `T` using `T=cell(mult^2,1)`.
- Lines 48: computes `T((n^2+1):((n+1)^2))` using `T((n^2+1):((n+1)^2))=irr_sph_ten(mult,n)`.
- Lines 65: computes `L` using `L=pauli(mult)`.
- Lines 71: computes `T{1}` using `T{1}=((-1)^k)*(2^(-k/2))*L.p^k`.
- Lines 75: computes `q` using `q=k-n+2; T{n}=(1/sqrt((k+q)*(k-q+1)))*(L.m*T{n-1}-T{n-1}*L.m)`.

### Local helper functions

- Line 90: `grumble()` — `function grumble(mult,k)`.
  - Representative operation: `if (~isnumeric(mult))||(~isscalar(mult))||(~isreal(mult))|| (~isfinite(mult))||(mult<1)||(mod(mult,1)~=0)`.
  - Representative operation: `(~isfinite(mult))||(mult<1)||(mod(mult,1)~=0)`.

## Syntax

```matlab
T=irr_sph_ten(mult,k)
```

## Parameters / inputs

- mult -multiplicity of the spin in question
- k -irreducible spherical tensor rank (optional)

## Outputs

- T -a two-argument call returns a cell array of
- tensors of rank k in the order of decreasing
- projection. A single argument call produces
- tensors of all ranks and puts them into a
- cell array in the order of increasing rank,
- and decreasing projection within each rank.
- Note: operator normalisation in spin dynamics is not a good
- idea. The only way to make the formalism independent of the
- spin quantum number is to impose identical commutation rela-
- tions rather than equal matrix norms.

## Implementation structure

- Single-spin irreducible spherical tensor operators T(k,m)
- obeying the following commutation relation:
- [Lz,T_km]=m*T_km
- T=irr_sph_ten(mult,k)
- mult -multiplicity of the spin in question
- k -irreducible spherical tensor rank (optional)
- T -a two-argument call returns a cell array of
- tensors of rank k in the order of decreasing
- projection. A single argument call produces
- tensors of all ranks and puts them into a
- cell array in the order of increasing rank,
- and decreasing projection within each rank.

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `grumble()`, `speye()`, `pauli()`, `isscalar()`.
