# kernel/overloads/@ttclass/ranks.m

- Signature: `ttranks=ranks(ttrain)`

## Purpose

Returns the bond dimensions of a tensor train. Syntax: ttranks=ranks(ttrain)

## Physical / mathematical content

- Tensor-train linear algebra. These files implement compressed high-dimensional operators and AMEn/SVD-based algebra in tensor-train format.
- The optimisation logic is Newton or Newton-like: search directions use first- and second-order local curvature information, usually with regularisation or line-search safeguards.

## Numerical / algorithmic content

## Parameters / inputs

- ttrain -a tensor train object

## Outputs

- ttranks -(ncores+1) by (ntrains) array; the
- first and the last elements for each
- train are 1

## Implementation structure

- Returns the bond dimensions of a tensor train. Syntax:
- ttranks=ranks(ttrain)
- ttrain -a tensor train object
- ttranks -(ncores+1) by (ntrains) array; the
- first and the last elements for each
- train are 1
- Get core array dimensions
- Preallocate the answer
- Loop over the buffer
- Extract the ranks
- I refrain from publishing for fear that disputes and controversies
- may be raised against me by ignoramuses.
