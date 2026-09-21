# kernel/overloads/@ttclass/size.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/kernel/overloads/@ttclass/size.m`
- Signature: `varargout=size(tt,dim)`
- Total lines: 61

## Purpose

Returns the size of the matrix represented by a tensor train. The output mimics Matlab's size function. Syntax: [m,n]=size(tt,dim)

## Physical / mathematical content

- Tensor-train linear algebra. These files implement compressed high-dimensional operators and AMEn/SVD-based algebra in tensor-train format.

## Numerical / algorithmic content

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
