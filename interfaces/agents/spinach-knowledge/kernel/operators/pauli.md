# kernel/operators/pauli.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/kernel/operators/pauli.m`
- Signature: `S=pauli(mult)`
- Total lines: 98

## Purpose

Pauli spin operators (sparse, see below for normalisa- tion conventions) for a spin of a user-specified ener- gy level multiplicity. Syntax: S=pauli(mult)

## Physical / mathematical content

- Operator-construction utilities. They build bases and irreducible tensor representations for spin, bosonic, and transition operators.

## Numerical / algorithmic content

- The file contains an explicit `grumble(...)` validator, which is Spinach convention for front-loading dimension, type, and regime checks before expensive linear-algebra work begins.
- The file also defines local helper function(s): `grumble()`. This usually means the public entry point is supported by tightly coupled validation or helper logic kept private to the file.

## Code-derived implementation details

### Comment-guided execution stages

- Lines 47-48: Ensure internal consistency; implemented by `grumble(mult); mult=double(mult)`.
- Lines 50-51: Get Pauli matrices; implemented by `if mult==2`.
- Lines 53-54: Spin-half matrices are hard-coded for speed; implemented by `S.u=sparse(complex([1 0; 0 1]))`.
- Lines 61-62: Spin-one matrices are hard-coded for speed; implemented by `S.u=sparse(complex([1 0 0; 0 1 0; 0 0 1]))`.
- Lines 69-70: Everything else gets generated; implemented by `spin=(mult-1)/2; prjs=((mult-1):-1:0)-spin`.
- Lines 79-80: X and Y operators are combinations; implemented by `S.x=(S.p+S.m)/2; S.y=(S.p-S.m)/2i`.

### Control flow inferred from the code

- Line 51: conditional branch on `mult==2`.

### Key state/data transformations

- Lines 48: computes `grumble(mult); mult` using `grumble(mult); mult=double(mult)`.
- Lines 54: computes `S.u` using `S.u=sparse(complex([1 0; 0 1]))`.
- Lines 55: computes `S.p` using `S.p=sparse(complex([0 1; 0 0]))`.
- Lines 56: computes `S.m` using `S.m=sparse(complex([0 0; 1 0]))`.
- Lines 57: computes `S.z` using `S.z=sparse(complex([0.5 0; 0 -0.5]))`.
- Lines 70: computes `spin` using `spin=(mult-1)/2; prjs=((mult-1):-1:0)-spin`.
- Lines 80: computes `S.x` using `S.x=(S.p+S.m)/2; S.y=(S.p-S.m)/2i`.

### Local helper functions

- Line 85: `grumble()` — `function grumble(mult)`. Any refusal to recognize reality, for any reason whatever, has dis- astrous consequences. There are no evil thoughts except one --the
  - Representative operation: `if (~isnumeric(mult))||(~isreal(mult))|| (~isscalar(mult))||(mod(mult,1)~=0)||(mult<1)`.
  - Representative operation: `(~isscalar(mult))||(mod(mult,1)~=0)||(mult<1)`.

## Parameters / inputs

- mult -an integer specifying the
- multiplicity of the spin

## Outputs

- S.u -unit operator
- S.p -raising operator
- S.m -lowering operator
- S.x -Sx observable operator
- S.y -Sy observable operator
- S.z -Sz observable operator
- Note: the matrices are normalised to obey the following
- commutation relations for all multiplicities:
- [S.x,S.y]=1i*S.z
- [S.y,S.z]=1i*S.x
- [S.z,S.x]=1i*S.y
- Note: raising and lowering operators are defined as:
- S.p=S.x+1i*S.y
- S.m=S.x-1i*S.y
- Note: arrays are declared complex at creation to avoid
- expensive reallocation operations later on.

## Implementation structure

- Pauli spin operators (sparse, see below for normalisa-
- tion conventions) for a spin of a user-specified ener-
- gy level multiplicity. Syntax:
- S=pauli(mult)
- mult -an integer specifying the
- multiplicity of the spin
- S.u -unit operator
- S.p -raising operator
- S.m -lowering operator
- S.x -Sx observable operator
- S.y -Sy observable operator
- S.z -Sz observable operator

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `grumble()`, `double()`, `complex()`, `spdiags()`, `speye()`, `isscalar()`.
