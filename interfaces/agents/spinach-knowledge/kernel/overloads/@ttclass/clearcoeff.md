# kernel/overloads/@ttclass/clearcoeff.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/kernel/overloads/@ttclass/clearcoeff.m`
- Signature: `tt=clearcoeff(tt)`
- Total lines: 48

## Purpose

Absorbs physical coefficients into tensor train cores without changing the value represented by the tensor train. Syntax: tt=clearcoeff(tt)

## Physical / mathematical content

- Tensor-train linear algebra. These files implement compressed high-dimensional operators and AMEn/SVD-based algebra in tensor-train format.

## Numerical / algorithmic content


## Code-derived implementation details

### Comment-guided execution stages

- Lines 22-23: Get the number of cores and trains; implemented by `[ncores,ntrains]=size(tt.cores)`.
- Lines 25-26: Loop over the trains in the buffer; implemented by `for n=1:ntrains`.
- Lines 28-29: Scale the coefficient; implemented by `A=tt.coeff(1,n)^(1/ncores)`.
- Lines 31-32: Apply it to cores; implemented by `for k=1:ncores`.
- Lines 36-37: Erase the coefficient; implemented by `tt.coeff(1,n)=1`.

### Control flow inferred from the code

- Line 26: `for` loop over `n=1:ntrains`.
- Line 32: `for` loop over `k=1:ncores`.

### Key state/data transformations

- Lines 23: computes `[ncores,ntrains]` using `[ncores,ntrains]=size(tt.cores)`.
- Lines 29: computes `A` using `A=tt.coeff(1,n)^(1/ncores)`.
- Lines 33: computes `tt.cores{k,n}` using `tt.cores{k,n}=A*tt.cores{k,n}`.
- Lines 37: computes `tt.coeff(1,n)` using `tt.coeff(1,n)=1`.

## Parameters / inputs

- tt -tensor train object

## Outputs

- tt -tensor train object with each coefficient distributed
- into its cores and the coefficient array set to one

## Implementation structure

- Absorbs physical coefficients into tensor train cores without
- changing the value represented by the tensor train. Syntax:
- tt=clearcoeff(tt)
- tt -tensor train object
- tt -tensor train object with each coefficient distributed
- into its cores and the coefficient array set to one
- Get the number of cores and trains
- Loop over the trains in the buffer
- Scale the coefficient
- Apply it to cores
- Erase the coefficient
- "Moral outrage is a middle-class luxury."
