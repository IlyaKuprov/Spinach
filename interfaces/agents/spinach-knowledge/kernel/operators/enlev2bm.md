# kernel/operators/enlev2bm.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/kernel/operators/enlev2bm.m`
- Signature: `[states,coeffs]=enlev2bm(nlevels,lvl_num)`
- Total lines: 59

## Purpose

Bosonic monomial expansion of specific bosonic energy level projectors. Syntax: [states,coeffs]=enlev2bm(nlevels,lvl_num)

## Physical / mathematical content

- Operator-construction utilities. They build bases and irreducible tensor representations for spin, bosonic, and transition operators.

## Numerical / algorithmic content

- The file contains an explicit `grumble(...)` validator, which is Spinach convention for front-loading dimension, type, and regime checks before expensive linear-algebra work begins.
- The file also defines local helper function(s): `grumble()`. This usually means the public entry point is supported by tightly coupled validation or helper logic kept private to the file.

## Code-derived implementation details

### Comment-guided execution stages

- Lines 30-31: Check consistency; implemented by `grumble(nlevels,lvl_num)`.
- Lines 33-34: Energy level projector; implemented by `P=zeros(nlevels,nlevels)`.
- Lines 37-38: Bosonic monomial expansion; implemented by `[states,coeffs]=oper2bm(P)`.

### Key state/data transformations

- Lines 34: computes `P` using `P=zeros(nlevels,nlevels)`.
- Lines 35: computes `P(lvl_num,lvl_num)` using `P(lvl_num,lvl_num)=1`.
- Lines 38: computes `[states,coeffs]` using `[states,coeffs]=oper2bm(P)`.

### Local helper functions

- Line 43: `grumble()` — `function grumble(nlevels,lvl_num)`. Истребление зануд -долг каждого порядочного
  - Representative operation: `if (~isnumeric(nlevels))||(~isscalar(nlevels))|| (~isreal(nlevels))||(nlevels<1)`.
  - Representative operation: `(~isreal(nlevels))||(nlevels<1)`.

## Parameters / inputs

- nlevels -number of energy levels in the mode,
- a positive integer
- lvl_num -energy level number, counting from the
- empty mode state upwards

## Outputs

- states -states, in the Spinach BM basis index-
- ing, that contribute to the operator in
- question; use lin2kq to convert to K,Q
- bosonic monomial indices
- coeffs -coefficients with which the BMs enter
- the linear combination

## Implementation structure

- Bosonic monomial expansion of specific bosonic energy
- level projectors. Syntax:
- [states,coeffs]=enlev2bm(nlevels,lvl_num)
- nlevels -number of energy levels in the mode,
- a positive integer
- lvl_num -energy level number, counting from the
- empty mode state upwards
- states -states, in the Spinach BM basis index-
- ing, that contribute to the operator in
- question; use lin2kq to convert to K,Q
- bosonic monomial indices
- coeffs -coefficients with which the BMs enter

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `grumble()`, `oper2bm()`, `isscalar()`.
