# kernel/overloads/@ttclass/sizes.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/kernel/overloads/@ttclass/sizes.m`
- Signature: `modesizes=sizes(tt)`
- Total lines: 43

## Purpose

Returns mode sizes (physical dimensions of each core) of a tensor train. Syntax: modesizes=sizes(tt)

## Physical / mathematical content

- Tensor-train linear algebra. These files implement compressed high-dimensional operators and AMEn/SVD-based algebra in tensor-train format.

## Numerical / algorithmic content

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
