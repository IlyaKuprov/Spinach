# kernel/overloads/@polyadic/nnz.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/kernel/overloads/@polyadic/nnz.m`
- Signature: `answer=nnz(p)`
- Total lines: 45

## Purpose

Number of non-zeroes in all kernels of the polyadic. Syntax: answer=nnz(p)

## Physical / mathematical content

- Polyadic tensor-product linear algebra. The emphasis is compressed operator representation, deferred algebra, and efficient Kronecker-structured manipulations.

## Numerical / algorithmic content


## Code-derived implementation details

### Comment-guided execution stages

- Lines 19-20: Start from zero; implemented by `answer=0`.
- Lines 22-23: Loop over cores; implemented by `for n=1:numel(p.cores)`.
- Lines 29-30: Loop over prefix; implemented by `for n=1:numel(p.prefix)`.
- Lines 34-35: Loop over suffix; implemented by `for n=1:numel(p.suffix)`.

### Control flow inferred from the code

- Line 23: `for` loop over `n=1:numel(p.cores)`.
- Line 24: `for` loop over `k=1:numel(p.cores{n})`.
- Line 30: `for` loop over `n=1:numel(p.prefix)`.
- Line 35: `for` loop over `n=1:numel(p.suffix)`.

### Key state/data transformations

- Lines 20: computes `answer` using `answer=0`.

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
