# kernel/overloads/@ttclass/mldivide.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/kernel/overloads/@ttclass/mldivide.m`
- Signature: `x=mldivide(A,y)`
- Total lines: 50

## Purpose

Solves a linear system with tensor train objects. Syntax: x=mldivide(A,y)

## Physical / mathematical content

- Tensor-train linear algebra. These files implement compressed high-dimensional operators and AMEn/SVD-based algebra in tensor-train format.

## Numerical / algorithmic content

## Parameters / inputs

- A -ttclass matrix
- y -ttclass vector

## Outputs

- x -ttclass vector
- Note: the AMEn-solve algorithm is applied to symmetrised
- system (A'*A)*x=A'*y

## Implementation structure

- Solves a linear system with tensor train objects. Syntax:
- x=mldivide(A,y)
- A -ttclass matrix
- y -ttclass vector
- x -ttclass vector
- Note: the AMEn-solve algorithm is applied to symmetrised
- system (A'*A)*x=A'*y
- Shrink the operands
- Form a symmetrised system
- Solve it with AMEn algorithm
- Complain and bomb out
- The penalty for success is to be bored by the people

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `shrink()`, `amensolve()`.
