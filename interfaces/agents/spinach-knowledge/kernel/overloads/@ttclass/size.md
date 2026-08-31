# kernel/overloads/@ttclass/size.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/kernel/overloads/@ttclass/size.m`
- Signature: `varargout=size(tt,dim)`
- Total lines: 61

## Purpose

Returns the size of the matrix represented by a tensor train. The output mimics Matlab's size function. Syntax: [m,n]=size(tt,dim)

## Physical / mathematical content

- Tensor-train linear algebra. These files implement compressed high-dimensional operators and AMEn/SVD-based algebra in tensor-train format.

## Numerical / algorithmic content


## Code-derived implementation details

### Comment-guided execution stages

- Lines 26-27: Multiply up physical dimensions of all cores; implemented by `m=prod(cellfun(@(x)size(x,2),tt.cores),1)`.
- Lines 30-31: Compose the answer; implemented by `if (nargin==1)&&(nargout<=1)`.
- Lines 44-45: Check for infinities; implemented by `if any(cellfun(@(x)any(x(:)>intmax),varargout))`.

### Control flow inferred from the code

- Line 31: conditional branch on `(nargin==1)&&(nargout<=1)`.
- Line 45: conditional branch on `any(cellfun(@(x)any(x(:)>intmax),varargout))`.

### Key state/data transformations

- Lines 27: computes `m` using `m=prod(cellfun(@(x)size(x,2),tt.cores),1)`.
- Lines 28: computes `n` using `n=prod(cellfun(@(x)size(x,3),tt.cores),1)`.
- Lines 32: computes `varargout{1}` using `varargout{1}=[m(1) n(1)]`.
- Lines 35: computes `varargout{2}` using `varargout{2}=n(1)`.

## Parameters / inputs

- tt -a tensor train representation of a matrix
- dim -(optional) integer specifying dimension

## Outputs

- m,n -integers specifying the dimensions
- Note: for large tt matrices, m and n may be too large to fit
- into the maximum integer permitted by Matlab.

## Implementation structure

- Returns the size of the matrix represented by a tensor train.
- The output mimics Matlab's size function. Syntax:
- [m,n]=size(tt,dim)
- tt -a tensor train representation of a matrix
- dim -(optional) integer specifying dimension
- m,n -integers specifying the dimensions
- Note: for large tt matrices, m and n may be too large to fit
- into the maximum integer permitted by Matlab.
- Multiply up physical dimensions of all cores
- Compose the answer
- Check for infinities
- Will fluorine ever have practical applications? It is very

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `cellfun()`, `elseif()`, `any()`.
