# kernel/operators/bos2ist.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/kernel/operators/bos2ist.m`
- Signature: `[states,coeffs]=bos2ist(prod_spec,nlevels)`
- Total lines: 86

## Purpose

Irreducible spherical tensor expansion of a user-specified bosonic operator product. Syntax: [states,coeffs]=bos2ist(prod_spec,lvl_num)

## Physical / mathematical content

- Operator-construction utilities. They build bases and irreducible tensor representations for spin, bosonic, and transition operators.

## Numerical / algorithmic content

- The file contains an explicit `grumble(...)` validator, which is Spinach convention for front-loading dimension, type, and regime checks before expensive linear-algebra work begins.
- The file also defines local helper function(s): `grumble()`. This usually means the public entry point is supported by tightly coupled validation or helper logic kept private to the file.

## Code-derived implementation details

### Comment-guided execution stages

- Lines 32-33: Check consistency; implemented by `grumble(prod_spec,nlevels)`.
- Lines 35-36: Weyl operators; implemented by `W=weyl(nlevels)`.
- Lines 38-39: Compute the product; implemented by `P=speye(nlevels,nlevels)`.
- Lines 45-46: Creation; implemented by `P=P*W.c`.
- Lines 50-51: Annihilation; implemented by `P=P*W.a`.
- Lines 55-56: Number; implemented by `P=P*W.n`.
- Lines 60-61: Complain and bomb out; implemented by `error('unknown operator specification.')`.
- Lines 66-67: Spherical tensor expansion; implemented by `[states,coeffs]=oper2ist(P)`.

### Control flow inferred from the code

- Line 40: `for` loop over `n=1:numel(prod_spec)`.
- Line 41: dispatches on `prod_spec(n)`; cases `'C'`, `'A'`, `'N'`.

### Key state/data transformations

- Lines 36: computes `W` using `W=weyl(nlevels)`.
- Lines 39: computes `P` using `P=speye(nlevels,nlevels)`.
- Lines 67: computes `[states,coeffs]` using `[states,coeffs]=oper2ist(P)`.

### Local helper functions

- Line 72: `grumble()` — `function grumble(prod_spec,nlevels)`. Wenn sich der Most auch ganz absurd gebärdet, Es gibt zuletzt doch noch n’ Wein.
  - Representative operation: `if ~char(prod_spec)`.
  - Representative operation: `error('prod_spec must be a character string.')`.

## Parameters / inputs

- prod_spec -bosonic operator product specification
- in which 'C' stands for creation opera-
- tor and 'A' for annihilation operator,
- for example 'CCAA'
- nlevels -number of energy levels in the trunca-
- ted bosonic mode

## Outputs

- states -states, in the Spinach IST basis index-
- ing, that contribute to the operator in
- question; use lin2lm to convert to L,M
- spherical tensor indices
- coeffs -coefficients with which the ISTs enter
- the linear combination

## Implementation structure

- Irreducible spherical tensor expansion of a user-specified
- bosonic operator product. Syntax:
- [states,coeffs]=bos2ist(prod_spec,lvl_num)
- prod_spec -bosonic operator product specification
- in which 'C' stands for creation opera-
- tor and 'A' for annihilation operator,
- for example 'CCAA'
- nlevels -number of energy levels in the trunca-
- ted bosonic mode
- states -states, in the Spinach IST basis index-
- ing, that contribute to the operator in
- question; use lin2lm to convert to L,M

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `grumble()`, `weyl()`, `speye()`, `prod_spec()`, `oper2ist()`, `char()`, `isscalar()`.
