# kernel/overloads/@ttclass/sum.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/kernel/overloads/@ttclass/sum.m`
- Signature: `answer=sum(ttrain,dim)`
- Total lines: 82

## Purpose

Sum of elements of a tensor train representation of a matrix. Syntax: answer=sum(ttrain,dim)

## Physical / mathematical content

- Tensor-train linear algebra. These files implement compressed high-dimensional operators and AMEn/SVD-based algebra in tensor-train format.
- The relevant state manifold is the singlet/triplet decomposition, where permutation symmetry controls selection rules, relaxation susceptibility, and convertibility to ordinary magnetisation.

## Numerical / algorithmic content


## Code-derived implementation details

### Comment-guided execution stages

- Lines 25-26: Get sizes and ranks; implemented by `[ncores,ntrains]=size(ttrain.cores)`.
- Lines 30-31: If all dimensions are singleton, return a scalar immediately; implemented by `if all(tt_sizes(:)==1), answer=full(ttrain); return; end`.
- Lines 33-35: In dim is omitted, choose first non-singleton dimension (this mimics the Matlab behaviour for matices); implemented by `if nargin==1`.
- Lines 44-45: Make an auxiliary tensor train; implemented by `answer=ttclass`.
- Lines 50-51: Run through all tensor trains; implemented by `for n=1:ntrains`.
- Lines 53-54: Run through all cores; implemented by `for k=1:ncores`.
- Lines 56-57: Sum the appropriate dimension in each core; implemented by `answer.cores{k,n}=sum(ttrain.cores{k,n},dim+1)`.
- Lines 59-60: Reshape the core; implemented by `switch dim`.
- Lines 73-74: If all modes are singleton, transform to a scalar; implemented by `if all(all(sizes(answer)==1)), answer=full(answer); end`.

### Control flow inferred from the code

- Line 31: conditional branch on `all(tt_sizes(:)==1), answer=full(ttrain); return; end`.
- Line 35: conditional branch on `nargin==1`.
- Line 37: `for` loop over `k=2:-1:1`.
- Line 38: conditional branch on `~all(tt_sizes(:,k)==1)`.
- Line 51: `for` loop over `n=1:ntrains`.
- Line 54: `for` loop over `k=1:ncores`.
- Line 60: dispatches on `dim`; cases `1`, `2`.
- Line 74: conditional branch on `all(all(sizes(answer)==1)), answer=full(answer); end`.

### Key state/data transformations

- Lines 26: computes `[ncores,ntrains]` using `[ncores,ntrains]=size(ttrain.cores)`.
- Lines 27: computes `tt_ranks` using `tt_ranks=ranks(ttrain)`.
- Lines 28: computes `tt_sizes` using `tt_sizes=sizes(ttrain)`.
- Lines 36: computes `dim` using `dim=0`.
- Lines 45: computes `answer` using `answer=ttclass`.
- Lines 46: computes `answer.coeff` using `answer.coeff=ttrain.coeff`.
- Lines 47: computes `answer.tolerance` using `answer.tolerance=zeros(1,ntrains)`.
- Lines 48: computes `answer.cores` using `answer.cores=cell(ncores,ntrains)`.
- Lines 57: computes `answer.cores{k,n}` using `answer.cores{k,n}=sum(ttrain.cores{k,n},dim+1)`.

## Parameters / inputs

- ttrain -tensor train object representing
- a matrix
- dim -summation dimension, 1 or 2

## Outputs

- answer -tensor train or flat representation
- of the summation result

## Implementation structure

- Sum of elements of a tensor train representation of
- a matrix. Syntax:
- answer=sum(ttrain,dim)
- ttrain -tensor train object representing
- a matrix
- dim -summation dimension, 1 or 2
- answer -tensor train or flat representation
- of the summation result
- Get sizes and ranks
- If all dimensions are singleton, return a scalar immediately
- In dim is omitted, choose first non-singleton dimension
- (this mimics the Matlab behaviour for matices)

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `ranks()`, `sizes()`, `all()`, `tt_sizes()`, `tt_ranks()`.
