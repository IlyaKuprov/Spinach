# kernel/overloads/@ttclass/rdivide.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/kernel/overloads/@ttclass/rdivide.m`
- Signature: `a=rdivide(a,b)`
- Total lines: 42

## Purpose

Divides a tensor train object by a scalar. Syntax: c=rdivide(a,b)

## Physical / mathematical content

- Tensor-train linear algebra. These files implement compressed high-dimensional operators and AMEn/SVD-based algebra in tensor-train format.

## Numerical / algorithmic content


## Code-derived implementation details

### Comment-guided execution stages

- Lines 22-23: Division of tensor train by a scalar; implemented by `if isa(a,'ttclass')&&isscalar(b)`.
- Lines 25-26: Divide the coefficients and update the tolerances; implemented by `a.coeff=a.coeff/b; a.tolerance=a.tolerance/abs(b)`.
- Lines 30-31: Complain and bomb out; implemented by `error('the first argument must be a tensor train and the second one a scalar.')`.

### Control flow inferred from the code

- Line 23: conditional branch on `isa(a,'ttclass')&&isscalar(b)`.

### Key state/data transformations

- Lines 26: computes `a.coeff` using `a.coeff=a.coeff/b; a.tolerance=a.tolerance/abs(b)`.

## Parameters / inputs

- a -a ttclass object
- b -a numeric scalar

## Outputs

- c -a ttclass object

## Implementation structure

- Divides a tensor train object by a scalar. Syntax:
- c=rdivide(a,b)
- a -a ttclass object
- b -a numeric scalar
- c -a ttclass object
- Division of tensor train by a scalar
- Divide the coefficients and update the tolerances
- Complain and bomb out
- Documentation is like sex: when it is good, it is
- very, very good, and when it is bad it's still bet-
- ter than nothing.
- Jim Hargrove

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `isscalar()`.
