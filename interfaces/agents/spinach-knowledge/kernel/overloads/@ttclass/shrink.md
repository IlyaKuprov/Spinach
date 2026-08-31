# kernel/overloads/@ttclass/shrink.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/kernel/overloads/@ttclass/shrink.m`
- Signature: `ttrain=shrink(ttrain)`
- Total lines: 48

## Purpose

Approximates a given tensor train with lower TT-ranks. Syntax: ttrain=shrink(ttrain)

## Physical / mathematical content

- Tensor-train linear algebra. These files implement compressed high-dimensional operators and AMEn/SVD-based algebra in tensor-train format.

## Numerical / algorithmic content


## Code-derived implementation details

### Comment-guided execution stages

- Lines 21-22: Read train sizes; implemented by `[d,~]=size(ttrain.cores)`.
- Lines 24-25: Summation; implemented by `ttrain=pack(ttrain)`.
- Lines 27-28: Left-to-right orthogonalisation; implemented by `ttrain=ttort(ttrain,+1)`.
- Lines 30-31: Check the norm and escape if the object is zero; implemented by `nrm=ttrain.coeff*norm(ttrain.cores{d,1}(:),2)`.
- Lines 34-35: Truncation; implemented by `ttrain=truncate(ttrain)`.
- Lines 37-39: Convert to a scalar if appropriate; implemented by `if all(all(cellfun(@(x)size(x,2),ttrain.cores)==1))&& all(all(cellfun(@(x)size(x,3),ttrain.cores)==1))`.

### Control flow inferred from the code

- Line 32: conditional branch on `nrm==0, ttrain=0*unit_like(ttrain); return; end`.
- Line 38: conditional branch on `all(all(cellfun(@(x)size(x,2),ttrain.cores)==1))&&`.

### Key state/data transformations

- Lines 22: computes `[d,~]` using `[d,~]=size(ttrain.cores)`.
- Lines 25: computes `ttrain` using `ttrain=pack(ttrain)`.
- Lines 31: computes `nrm` using `nrm=ttrain.coeff*norm(ttrain.cores{d,1}(:),2)`.

## Parameters / inputs

- ttrain -a tensor train object

## Outputs

- ttrain -compressed tensor train object with
- right-to-left orthogonalisation

## Implementation structure

- Approximates a given tensor train with lower TT-ranks. Syntax:
- ttrain=shrink(ttrain)
- ttrain -a tensor train object
- ttrain -compressed tensor train object with
- right-to-left orthogonalisation
- Read train sizes
- Summation
- Left-to-right orthogonalisation
- Check the norm and escape if the object is zero
- Truncation
- Convert to a scalar if appropriate
- It's a tough life, being small and delicious.

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `pack()`, `ttort()`, `unit_like()`, `truncate()`, `all()`, `cellfun()`.
