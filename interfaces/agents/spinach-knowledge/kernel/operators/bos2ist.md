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
