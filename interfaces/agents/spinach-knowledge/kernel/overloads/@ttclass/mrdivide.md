# kernel/overloads/@ttclass/mrdivide.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/kernel/overloads/@ttclass/mrdivide.m`
- Signature: `a=mrdivide(a,b)`
- Total lines: 41

## Purpose

Divides a tensor train object by a scalar. Syntax: c=mrdivide(a,b)

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

- a -tensor train object
- b -a scalar

## Outputs

- c -tensor train object

## Implementation structure

- Divides a tensor train object by a scalar. Syntax:
- c=mrdivide(a,b)
- a -tensor train object
- b -a scalar
- c -tensor train object
- Division of tensor train by a scalar
- Divide the coefficients and update the tolerances
- Complain and bomb out
- It is dangerous to be right in matters on which the established
- authorities are wrong.
- Voltaire

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `isscalar()`.
