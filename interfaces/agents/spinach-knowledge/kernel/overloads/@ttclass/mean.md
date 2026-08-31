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


## Code-derived implementation details

### Comment-guided execution stages

- Lines 23-24: Get sizes and ranks; implemented by `[ncores,ntrains]=size(ttrain.cores)`.
- Lines 28-29: If all dimensions are singleton, return a scalar immediately; implemented by `if all(tt_sizes(:)==1), answer=full(ttrain); return; end`.
- Lines 31-33: In dim is omitted, choose first non-singleton dimension (this mimics the Matlab behaviour for matices); implemented by `if nargin==1`.
- Lines 42-43: Make an auxiliary tensor train; implemented by `answer=ttclass`.
- Lines 48-49: Run through all tensor trains; implemented by `for n=1:ntrains`.
- Lines 51-52: Run through all cores; implemented by `for k=1:ncores`.
- Lines 54-55: Sum the appropriate dimension in each core; implemented by `answer.cores{k,n}=sum(ttrain.cores{k,n},dim+1)`.
- Lines 57-58: Reshape the core and divide the result by the correspondent dimension; implemented by `switch dim`.
- Lines 71-72: If all modes are singleton, transform to a scalar; implemented by `if all(all(sizes(answer)==1)), answer=full(answer); end`.

### Control flow inferred from the code

- Line 29: conditional branch on `all(tt_sizes(:)==1), answer=full(ttrain); return; end`.
- Line 33: conditional branch on `nargin==1`.
- Line 35: `for` loop over `k=2:-1:1`.
- Line 36: conditional branch on `~all(tt_sizes(:,k)==1)`.
- Line 49: `for` loop over `n=1:ntrains`.
- Line 52: `for` loop over `k=1:ncores`.
- Line 58: dispatches on `dim`; cases `1`, `2`.
- Line 72: conditional branch on `all(all(sizes(answer)==1)), answer=full(answer); end`.

### Key state/data transformations

- Lines 24: computes `[ncores,ntrains]` using `[ncores,ntrains]=size(ttrain.cores)`.
- Lines 25: computes `tt_ranks` using `tt_ranks=ranks(ttrain)`.
- Lines 26: computes `tt_sizes` using `tt_sizes=sizes(ttrain)`.
- Lines 34: computes `dim` using `dim=0`.
- Lines 43: computes `answer` using `answer=ttclass`.
- Lines 44: computes `answer.coeff` using `answer.coeff=ttrain.coeff`.
- Lines 45: computes `answer.tolerance` using `answer.tolerance=zeros(1,ntrains)`.
- Lines 46: computes `answer.cores` using `answer.cores=cell(ncores,ntrains)`.
- Lines 55: computes `answer.cores{k,n}` using `answer.cores{k,n}=sum(ttrain.cores{k,n},dim+1)`.

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
