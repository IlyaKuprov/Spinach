# kernel/overloads/@ttclass/sizes.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/kernel/overloads/@ttclass/sizes.m`
- Signature: `modesizes=sizes(tt)`
- Total lines: 43

## Purpose

Returns mode sizes (physical dimensions of each core) of a tensor train. Syntax: modesizes=sizes(tt)

## Physical / mathematical content

- Tensor-train linear algebra. These files implement compressed high-dimensional operators and AMEn/SVD-based algebra in tensor-train format.

## Numerical / algorithmic content


## Code-derived implementation details

### Comment-guided execution stages

- Lines 22-23: Determine the number of cores; implemented by `ncores=size(tt.cores,1)`.
- Lines 25-26: Preallocate the answer; implemented by `modesizes=zeros(ncores,2)`.
- Lines 28-29: Fill in the answer; implemented by `for k=1:ncores`.

### Control flow inferred from the code

- Line 29: `for` loop over `k=1:ncores`.

### Key state/data transformations

- Lines 23: computes `ncores` using `ncores=size(tt.cores,1)`.
- Lines 26: computes `modesizes` using `modesizes=zeros(ncores,2)`.
- Lines 30: computes `modesizes(k,1)` using `modesizes(k,1)=size(tt.cores{k,1},2)`.
- Lines 31: computes `modesizes(k,2)` using `modesizes(k,2)=size(tt.cores{k,1},3)`.

## Parameters / inputs

- tt -tensor train object

## Outputs

- modesizes -ncores by 2 array of physical dimensions
- of tensor train cores

## Implementation structure

- Returns mode sizes (physical dimensions of each core) of
- a tensor train. Syntax:
- modesizes=sizes(tt)
- tt -tensor train object
- modesizes -ncores by 2 array of physical dimensions
- of tensor train cores
- Determine the number of cores
- Preallocate the answer
- Fill in the answer
- Computer models are no different from fashion models: seductive,
- unreliable, easily corrupted, and they lead sensible people to
- make fools of themselves.

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `modesizes()`.
