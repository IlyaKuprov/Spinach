# kernel/overloads/@ttclass/isreal.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/kernel/overloads/@ttclass/isreal.m`
- Signature: `answer=isreal(tt)`
- Total lines: 51

## Purpose

Returns TRUE for real-valued tensor train objects. Syntax: answer=isreal(tt)

## Physical / mathematical content

- Tensor-train linear algebra. These files implement compressed high-dimensional operators and AMEn/SVD-based algebra in tensor-train format.

## Numerical / algorithmic content


## Code-derived implementation details

### Comment-guided execution stages

- Lines 20-21: Non-empty tensor trains should return true(); implemented by `if isa(tt,'ttclass')`.
- Lines 23-24: Check coefficient first; implemented by `answer=all(isreal(tt.coeff))`.
- Lines 26-27: If the coefficients are real, check the cores; implemented by `if answer`.
- Lines 38-39: Complain and bomb out; implemented by `error('input is not a ttclass.')`.

### Control flow inferred from the code

- Line 21: conditional branch on `isa(tt,'ttclass')`.
- Line 27: conditional branch on `answer`.
- Line 28: `for` loop over `n=1:tt.ntrains`.
- Line 29: `for` loop over `k=1:tt.ncores`.
- Line 31: conditional branch on `~answer, return; end`.

### Key state/data transformations

- Lines 24: computes `answer` using `answer=all(isreal(tt.coeff))`.

## Parameters / inputs

- tt -tensor train object

## Outputs

- answer -logical true when all coefficients and core
- elements of the tensor train are real

## Implementation structure

- Returns TRUE for real-valued tensor train objects. Syntax:
- answer=isreal(tt)
- tt -tensor train object
- answer -logical true when all coefficients and core
- elements of the tensor train are real
- Non-empty tensor trains should return true()
- Check coefficient first
- If the coefficients are real, check the cores
- Complain and bomb out
- Democracy is a pathetic belief in the collective wisdom
- of individual ignorance.
- H.L. Mencken

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `all()`.
