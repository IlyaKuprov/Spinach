# kernel/overloads/@ttclass/norm.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/kernel/overloads/@ttclass/norm.m`
- Signature: `ttnorm=norm(ttrain,norm_type) %#NORMOK`
- Total lines: 74

## Purpose

Computes the norm of the matrix represented by a tensor train. Syntax: ttnorm=norm(ttrain,norm_type)

## Physical / mathematical content

- Tensor-train linear algebra. These files implement compressed high-dimensional operators and AMEn/SVD-based algebra in tensor-train format.

## Numerical / algorithmic content

## Parameters / inputs

- ttrain -a tensor train representation of a matrix
- norm_type:
- norm_type=1 not available for ttclass
- norm_type=inf not available for ttclass
- norm_type=2 not available for ttclass
- norm_type='fro' returns the Frobenius norm

## Outputs

- ttnorm -a positive real number
- Note: only Frobenius norm is currently available for tensor trains;
- other norm types raise errors.

## Implementation structure

- Computes the norm of the matrix represented by a tensor train. Syntax:
- ttnorm=norm(ttrain,norm_type)
- ttrain -a tensor train representation of a matrix
- norm_type:
- norm_type=1 not available for ttclass
- norm_type=inf not available for ttclass
- norm_type=2 not available for ttclass
- norm_type='fro' returns the Frobenius norm
- ttnorm -a positive real number
- Note: only Frobenius norm is currently available for tensor trains;
- other norm types raise errors.
- Compute the norm

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `pack()`, `ttort()`.
