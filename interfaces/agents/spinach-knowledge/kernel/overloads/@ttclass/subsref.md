# kernel/overloads/@ttclass/subsref.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/kernel/overloads/@ttclass/subsref.m`
- Signature: `answer=subsref(ttrain,reference)`
- Total lines: 126

## Purpose

Dot and bracket property specifications for the tensor train class.

## Physical / mathematical content

- Tensor-train linear algebra. These files implement compressed high-dimensional operators and AMEn/SVD-based algebra in tensor-train format.

## Numerical / algorithmic content

- The file contains an explicit `grumble(...)` validator, which is Spinach convention for front-loading dimension, type, and regime checks before expensive linear-algebra work begins.
- The file also defines local helper function(s): `grumble()`, `ttclass_ind2sub()`. This usually means the public entry point is supported by tightly coupled validation or helper logic kept private to the file.

## Code-derived implementation details

### Comment-guided execution stages

- Lines 26-27: Methods and properties; implemented by `case '.'`.
- Lines 29-30: Return the output requested; implemented by `switch reference(1).subs`.
- Lines 41-42: Matrix element extraction; implemented by `case '()'`.
- Lines 44-45: Start with zero; implemented by `answer=0`.
- Lines 47-48: Convert indices; implemented by `if numel(reference(1).subs)~=2`.
- Lines 62-63: Multiply up the tensor train; implemented by `[ncores,ntrains]=size(ttrain.cores)`.
- Lines 74-75: Complain and bomb out; implemented by `error('unknown subscript reference type.')`.
- Lines 79-80: Allow nested indexing; implemented by `if numel(reference)>1`.

### Control flow inferred from the code

- Line 24: dispatches on `reference(1).type`; cases `'.'`, `'ncores', answer=ttrain.ncores`, `'ntrains', answer=ttrain.ntrains`, `'sizes', answer=ttrain.sizes`, `'ranks', answer=ttrain.ranks`, `'coeff', answer=ttrain.coeff`, `'cores', answer=ttrain.cores`, `'tolerance', answer=ttrain.tolerance`.
- Line 30: dispatches on `reference(1).subs`; cases `'ncores', answer=ttrain.ncores`, `'ntrains', answer=ttrain.ntrains`, `'sizes', answer=ttrain.sizes`, `'ranks', answer=ttrain.ranks`, `'coeff', answer=ttrain.coeff`, `'cores', answer=ttrain.cores`, `'tolerance', answer=ttrain.tolerance`.
- Line 48: conditional branch on `numel(reference(1).subs)~=2`.
- Line 53: conditional branch on `islogical(row_ind), row_ind=double(row_ind); end`.
- Line 54: conditional branch on `islogical(col_ind), col_ind=double(col_ind); end`.
- Line 64: `for` loop over `n=1:ntrains`.
- Line 66: `for` loop over `k=(ncores-1):(-1):1`.
- Line 80: conditional branch on `numel(reference)>1`.

### Key state/data transformations

- Lines 45: computes `answer` using `answer=0`.
- Lines 51: computes `siz` using `siz=sizes(ttrain)`.
- Lines 52: computes `row_ind` using `row_ind=reference(1).subs{1}; col_ind=reference(1).subs{2}`.
- Lines 56: computes `ind` using `ind=ttclass_ind2sub(siz(:,1),row_ind)`.
- Lines 57: computes `jnd` using `jnd=ttclass_ind2sub(siz(:,2),col_ind)`.
- Lines 63: computes `[ncores,ntrains]` using `[ncores,ntrains]=size(ttrain.cores)`.
- Lines 65: computes `x` using `x=ttrain.cores{ncores,n}(:,ind(ncores),jnd(ncores),:)`.

### Local helper functions

- Line 87: `grumble()` — `function grumble(row_ind,col_ind,siz)`.
  - Representative operation: `if (~isnumeric(row_ind))||(~isreal(row_ind))|| (mod(row_ind,1)~=0)||(row_ind<1)||(row_ind>flintmax)`.
  - Representative operation: `(mod(row_ind,1)~=0)||(row_ind<1)||(row_ind>flintmax)`.
- Line 113: `ttclass_ind2sub()` — `function ivec=ttclass_ind2sub(siz,ind)`. "A cactus is a very disappointed cucumber." A Russian saying
  - Representative operation: `d=numel(siz); ind=ind-1; ivec=zeros(d,1)`.
  - Representative operation: `for k=d:-1:1`.

## Syntax

```matlab
answer=subsref(ttrain,reference)
```

## Parameters / inputs

- ttrain -tensor train object
- reference -Matlab subscript reference structure

## Outputs

- answer -requested tensor train property, scalar
- matrix element, or nested subscript result

## Implementation structure

- Dot and bracket property specifications for the tensor train class.
- answer=subsref(ttrain,reference)
- ttrain -tensor train object
- reference -Matlab subscript reference structure
- answer -requested tensor train property, scalar
- matrix element, or nested subscript result
- Methods and properties
- Return the output requested
- Matrix element extraction
- Start with zero
- Convert indices
- Multiply up the tensor train

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `reference()`, `isscalar()`, `sizes()`, `islogical()`, `double()`, `grumble()`, `ttclass_ind2sub()`, `siz()`, `ind()`, `jnd()`, `ivec()`.
