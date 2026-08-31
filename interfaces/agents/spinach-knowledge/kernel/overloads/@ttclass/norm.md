# kernel/overloads/@ttclass/norm.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/kernel/overloads/@ttclass/norm.m`
- Signature: `ttnorm=norm(ttrain,norm_type) %#NORMOK`
- Total lines: 74

## Purpose

Computes the norm of the matrix represented by a tensor train. Syntax: ttnorm=norm(ttrain,norm_type)

## Physical / mathematical content

- Tensor-train linear algebra. These files implement compressed high-dimensional operators and AMEn/SVD-based algebra in tensor-train format.

## Numerical / algorithmic content


## Code-derived implementation details

### Comment-guided execution stages

- Lines 30-31: Compute the norm; implemented by `switch norm_type`.
- Lines 35-36: Frobenius norm; implemented by `ttrain=pack(ttrain); ttrain=ttort(ttrain,-1)`.
- Lines 41-42: Maximum absolute column sum; implemented by `error('1-norm is not available for ttclass')`.
- Lines 46-47: Maximum absolute row sum; implemented by `error('inf-norm is not available for ttclass')`.
- Lines 51-52: Maximum absolute eigenvalue; implemented by `error('2-norm is not available for ttclass')`.
- Lines 56-57: Complain and bomb out; implemented by `error('unrecognized norm type.')`.

### Control flow inferred from the code

- Line 31: dispatches on `norm_type`; cases `'fro'`, `1`, `inf`, `2`.

### Key state/data transformations

- Lines 36: computes `ttrain` using `ttrain=pack(ttrain); ttrain=ttort(ttrain,-1)`.
- Lines 37: computes `ttnorm` using `ttnorm=abs(ttrain.coeff)*norm(ttrain.cores{1,1}(:),2)`.

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
