# kernel/operators/centrans.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/kernel/operators/centrans.m`
- Signature: `A=centrans(mult,type)`
- Total lines: 87

## Purpose

Central transition operators of half-integer spins in the Pauli basis. Syntax: A=centrans(mult,type)

## Physical / mathematical content

- Operator-construction utilities. They build bases and irreducible tensor representations for spin, bosonic, and transition operators.

## Numerical / algorithmic content

- The file contains an explicit `grumble(...)` validator, which is Spinach convention for front-loading dimension, type, and regime checks before expensive linear-algebra work begins.
- The file also defines local helper function(s): `grumble()`. This usually means the public entry point is supported by tightly coupled validation or helper logic kept private to the file.

## Code-derived implementation details

### Comment-guided execution stages

- Lines 24-25: Check consistency; implemented by `grumble(mult,type)`.
- Lines 27-28: Build CT operator; implemented by `A=spalloc(mult,mult,2)`.
- Lines 33-34: Sx on central transition; implemented by `A(mult/2,mult/2+1)=0.5`.
- Lines 39-40: Sy on central transition; implemented by `A(mult/2,mult/2+1)=-0.5i`.
- Lines 45-46: Sz on central transition; implemented by `A(mult/2,mult/2)=0.5`.
- Lines 51-52: S+ on central transition; implemented by `A(mult/2,mult/2+1)=1`.
- Lines 56-57: S-on central transition; implemented by `A(mult/2+1,mult/2)=1`.
- Lines 61-62: Complain and bomb out; implemented by `error('unknown CT operator type.')`.
- Lines 66-67: Make complex; implemented by `A=complex(A)`.

### Control flow inferred from the code

- Line 29: dispatches on `type`; cases `'x'`, `'y'`, `'z'`, `'+'`, `'-'`.

### Key state/data transformations

- Lines 28: computes `A` using `A=spalloc(mult,mult,2)`.
- Lines 34: computes `A(mult/2,mult/2+1)` using `A(mult/2,mult/2+1)=0.5`.
- Lines 35: computes `A(mult/2+1,mult/2)` using `A(mult/2+1,mult/2)=0.5`.
- Lines 46: computes `A(mult/2,mult/2)` using `A(mult/2,mult/2)=0.5`.
- Lines 47: computes `A(mult/2+1,mult/2+1)` using `A(mult/2+1,mult/2+1)=-0.5`.

### Local helper functions

- Line 72: `grumble()` — `function grumble(mult,type)`. Prohibition agent Izzy Einstein bragged that he could find liquor in any city under 30 minutes. In Chicago it took him
  - Representative operation: `if (~isnumeric(mult))||(~isscalar(mult))|| (~isreal(mult))||(mult<2)||(mod(mult,2)~=0)`.
  - Representative operation: `(~isreal(mult))||(mult<2)||(mod(mult,2)~=0)`.

## Parameters / inputs

- mult -multipicity of the spin in question, an
- even positive integer
- type -operator type: 'z' for polarisation, '+'
- for raising, '-' for lowering

## Outputs

- A -central transition operator

## Implementation structure

- Central transition operators of half-integer spins in the
- Pauli basis. Syntax:
- A=centrans(mult,type)
- mult -multipicity of the spin in question, an
- even positive integer
- type -operator type: 'z' for polarisation, '+'
- for raising, '-' for lowering
- A -central transition operator
- Check consistency
- Build CT operator
- Sx on central transition
- Sy on central transition

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `grumble()`, `spalloc()`, `complex()`, `isscalar()`, `ischar()`, `ismember()`.
