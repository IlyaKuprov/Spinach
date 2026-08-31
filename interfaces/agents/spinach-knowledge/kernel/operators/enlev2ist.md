# kernel/operators/enlev2ist.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/kernel/operators/enlev2ist.m`
- Signature: `[states,coeffs]=enlev2ist(mult,lvl_num,particle)`
- Total lines: 92

## Purpose

Irreducible spherical tensor expansion of specific Zeeman energy level projectors. Syntax: [states,coeffs]=enlev2ist(mult,lvl_num,particle)

## Physical / mathematical content

- Operator-construction utilities. They build bases and irreducible tensor representations for spin, bosonic, and transition operators.

## Numerical / algorithmic content

- The file contains an explicit `grumble(...)` validator, which is Spinach convention for front-loading dimension, type, and regime checks before expensive linear-algebra work begins.
- The file also defines local helper function(s): `grumble()`. This usually means the public entry point is supported by tightly coupled validation or helper logic kept private to the file.

## Code-derived implementation details

### Comment-guided execution stages

- Lines 34-35: Check consistency; implemented by `grumble(mult,lvl_num)`.
- Lines 37-38: Energy level projector; implemented by `P=zeros(mult,mult)`.
- Lines 40-41: Energy level counting; implemented by `switch particle`.
- Lines 45-46: From the bottom up; implemented by `P(mult-lvl_num+1,mult-lvl_num+1)=1`.
- Lines 50-51: From top to bottom; implemented by `P(lvl_num,lvl_num)=1`.
- Lines 55-56: Complain and bomb out; implemented by `error('unknown particle type.')`.
- Lines 60-61: Spherical tensor expansion; implemented by `[states,coeffs]=oper2ist(P)`.

### Control flow inferred from the code

- Line 41: dispatches on `particle`; cases `'S'`, `'B'`.

### Key state/data transformations

- Lines 38: computes `P` using `P=zeros(mult,mult)`.
- Lines 46: computes `P(mult-lvl_num+1,mult-lvl_num+1)` using `P(mult-lvl_num+1,mult-lvl_num+1)=1`.
- Lines 51: computes `P(lvl_num,lvl_num)` using `P(lvl_num,lvl_num)=1`.
- Lines 61: computes `[states,coeffs]` using `[states,coeffs]=oper2ist(P)`.

### Local helper functions

- Line 66: `grumble()` — `function grumble(mult,lvl_num)`. |------------------------------------------------------------------------------------------------------------------|
  - Representative operation: `if (~isnumeric(mult))||(~isscalar(mult))|| (~isreal(mult))||(mult<1)`.
  - Representative operation: `(~isreal(mult))||(mult<1)`.

## Parameters / inputs

- mult -multipicity of the spin in question, a
- positive integer
- lvl_num -energy level number, counting from the
- bottom up for spins and from top down
- for bosons
- particle -particle type, 'S' for a spin and 'B'
- for a boson

## Outputs

- states -states, in the Spinach IST basis index-
- ing, that contribute to the operator in
- question; use lin2lm to convert to L,M
- spherical tensor indices
- coeffs -coefficients with which the ISTs enter
- the linear combination

## Implementation structure

- Irreducible spherical tensor expansion of specific Zeeman
- energy level projectors. Syntax:
- [states,coeffs]=enlev2ist(mult,lvl_num,particle)
- mult -multipicity of the spin in question, a
- positive integer
- lvl_num -energy level number, counting from the
- bottom up for spins and from top down
- for bosons
- particle -particle type, 'S' for a spin and 'B'
- for a boson
- states -states, in the Spinach IST basis index-
- ing, that contribute to the operator in

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `grumble()`, `oper2ist()`, `isscalar()`.
