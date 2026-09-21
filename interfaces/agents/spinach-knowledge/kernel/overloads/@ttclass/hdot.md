# kernel/overloads/@ttclass/hdot.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/kernel/overloads/@ttclass/hdot.m`
- Signature: `c=hdot(a,b)`
- Total lines: 71

## Purpose

Hadamard dot product between two tensor train matrices. Syntax: c=hdot(a,b)

## Physical / mathematical content

- Tensor-train linear algebra. These files implement compressed high-dimensional operators and AMEn/SVD-based algebra in tensor-train format.

## Numerical / algorithmic content

- The file contains an explicit `grumble(...)` validator, which is Spinach convention for front-loading dimension, type, and regime checks before expensive linear-algebra work begins.
- The file also defines local helper function(s): `size()`. This usually means the public entry point is supported by tightly coupled validation or helper logic kept private to the file.

## Parameters / inputs

- a,b -tensor train objects representing numerical
- arrays of the same dimensions and having
- the same internal topology

## Outputs

- c -Hadamard product of a and b, a scalar

## Implementation structure

- Hadamard dot product between two tensor train matrices. Syntax:
- c=hdot(a,b)
- a,b -tensor train objects representing numerical
- arrays of the same dimensions and having
- the same internal topology
- c -Hadamard product of a and b, a scalar
- Check consistency
- Read topology and initialize the answer
- Loop over TT buffers
- Multiply coefficients
- Loop over TT cores and compute dot product
- Add to the total

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `grumble()`, `ranks()`, `sizes()`, `conj()`, `ranks_b()`, `mode_sizes()`, `ranks_a()`, `all()`.
