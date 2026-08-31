# kernel/overloads/@rcv/minus.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/kernel/overloads/@rcv/minus.m`
- Signature: `A=minus(A,B)`
- Total lines: 51

## Purpose

Subtracts one RCV object from another. Syntax: A=minus(A,B)

## Physical / mathematical content

- RCV sparse-matrix storage utilities. The focus is data structure design for sparse linear algebra and low-overhead composition of large matrices.

## Numerical / algorithmic content

- The file contains an explicit `grumble(...)` validator, which is Spinach convention for front-loading dimension, type, and regime checks before expensive linear-algebra work begins.
- The file also defines local helper function(s): `grumble()`. This usually means the public entry point is supported by tightly coupled validation or helper logic kept private to the file.

## Code-derived implementation details

### Comment-guided execution stages

- Lines 21-22: Check consistency; implemented by `grumble(A,B)`.
- Lines 24-25: Just call plus; implemented by `A=plus(A,(-1)*B)`.

### Key state/data transformations

- Lines 25: computes `A` using `A=plus(A,(-1)*B)`.

### Local helper functions

- Line 30: `grumble()` — `function grumble(A,B)`. "My cat had been suffering from severe illness over the past month or so. This had meant that he had needed increasingly hands-on
  - Representative operation: `if (~isa(A,'rcv'))&&(~isa(B,'rcv'))`.
  - Representative operation: `error('at least one input must be an RCV sparse matrix.')`.

## Parameters / inputs

- A -left operand
- B -right operand

## Outputs

- A -result A-B as an RCV sparse matrix

## Implementation structure

- Subtracts one RCV object from another. Syntax:
- A=minus(A,B)
- A -left operand
- B -right operand
- A -result A-B as an RCV sparse matrix
- Check consistency
- Just call plus
- Consistency enforcement
- "My cat had been suffering from severe illness over the past month
- or so. This had meant that he had needed increasingly hands-on
- care. Due to a terminal diagnosis the decision to put him to
- sleep was made; that took place on Monday 11th April.

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `grumble()`, `plus()`.
