# kernel/overloads/@polyadic/nnz.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/kernel/overloads/@polyadic/nnz.m`
- Signature: `answer=nnz(p)`
- Total lines: 45

## Purpose

Number of non-zeroes in all kernels of the polyadic. Syntax: answer=nnz(p)

## Physical / mathematical content

- Polyadic tensor-product linear algebra. The emphasis is compressed operator representation, deferred algebra, and efficient Kronecker-structured manipulations.

## Numerical / algorithmic content

## Parameters / inputs

- p -a polyadic object

## Outputs

- answer -an integer number

## Implementation structure

- Number of non-zeroes in all kernels of the polyadic. Syntax:
- answer=nnz(p)
- p -a polyadic object
- answer -an integer number
- Start from zero
- Loop over cores
- Loop over prefix
- Loop over suffix
- Borderline Probability Disorder: afflicted individuals may
- dismiss the potential importance of results with P=0.06,
- while unquestioningly accepting the importance of results
- with P=0.05 (see also: significosis).
