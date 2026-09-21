# kernel/overloads/@ttclass/mean.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/kernel/overloads/@ttclass/mean.m`
- Signature: `answer=mean(ttrain,dim)`
- Total lines: 80

## Purpose

Mean of elements of a tensor train representation of a matrix. Syntax: answer=mean(ttrain,dim)

## Physical / mathematical content

- Tensor-train linear algebra. These files implement compressed high-dimensional operators and AMEn/SVD-based algebra in tensor-train format.
- The relevant state manifold is the singlet/triplet decomposition, where permutation symmetry controls selection rules, relaxation susceptibility, and convertibility to ordinary magnetisation.

## Numerical / algorithmic content

## Parameters / inputs

- ttrain -a tensor train representation of a matrix
- dim -dimension to operate on (dim=1 or dim=2)

## Outputs

- answer -the mean value computed along the speci-
- fied dimension

## Implementation structure

- Mean of elements of a tensor train representation of a
- matrix. Syntax:
- answer=mean(ttrain,dim)
- ttrain -a tensor train representation of a matrix
- dim -dimension to operate on (dim=1 or dim=2)
- answer -the mean value computed along the speci-
- fied dimension
- Get sizes and ranks
- If all dimensions are singleton, return a scalar immediately
- In dim is omitted, choose first non-singleton dimension
- (this mimics the Matlab behaviour for matices)
- Make an auxiliary tensor train

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `ranks()`, `sizes()`, `all()`, `tt_sizes()`, `tt_ranks()`.
