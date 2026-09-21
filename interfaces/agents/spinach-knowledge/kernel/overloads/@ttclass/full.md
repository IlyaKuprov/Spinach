# kernel/overloads/@ttclass/full.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/kernel/overloads/@ttclass/full.m`
- Signature: `answer=full(ttrain)`
- Total lines: 58

## Purpose

Converts a tensor train representation of a matrix into a matrix. Syntax: answer=full(ttrain)

## Physical / mathematical content

- Tensor-train linear algebra. These files implement compressed high-dimensional operators and AMEn/SVD-based algebra in tensor-train format.

## Numerical / algorithmic content

## Parameters / inputs

- ttrain -tensor train object

## Outputs

- answer -a full matrix
- Note: the result can be huge, careless use would crash the system.

## Implementation structure

- Converts a tensor train representation of a matrix
- into a matrix. Syntax:
- answer=full(ttrain)
- ttrain -tensor train object
- answer -a full matrix
- Note: the result can be huge, careless use would crash the system.
- Preallocate the result
- Get object dimensions
- Get tensor ranks
- Get mode sizes
- Loop over the buffer
- Multiply up the tensor train

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `ranks()`, `sizes()`, `ttranks()`, `modesizes()`.
