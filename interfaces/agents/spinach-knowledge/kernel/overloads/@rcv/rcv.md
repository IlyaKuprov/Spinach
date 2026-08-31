# kernel/overloads/@rcv/rcv.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/kernel/overloads/@rcv/rcv.m`
- Signature: `obj=rcv(varargin)`
- Total lines: 190

## Purpose

Creates an RCV (row-column-value storage) sparse matrix. Syntax: obj=rcv(M) obj=rcv(dim1,dim2) obj=rcv(R,C,V,dim1,dim2)

## Physical / mathematical content

- RCV sparse-matrix storage utilities. The focus is data structure design for sparse linear algebra and low-overhead composition of large matrices.

## Numerical / algorithmic content

- The implementation explicitly addresses performance engineering through parallel or GPU execution, which matters because Spinach operators can become extremely large after basis expansion or powder/spatial lifting.
- The file contains an explicit `grumble(...)` validator, which is Spinach convention for front-loading dimension, type, and regime checks before expensive linear-algebra work begins.
- The file also defines local helper function(s): `isnumeric()`, `ismatrix()`, `isfloat()`, `grumble()`. This usually means the public entry point is supported by tightly coupled validation or helper logic kept private to the file.

## Code-derived implementation details

### Comment-guided execution stages

- Lines 44-45: Check consistency; implemented by `grumble(varargin{:})`.
- Lines 49-50: Do nothing if already RCV; implemented by `if isa(varargin{1},'rcv')`.
- Lines 52-53: Return input; implemented by `obj=varargin{1}`.
- Lines 55-56: Convert Matlab matrices to RCV; implemented by `elseif ismatrix(varargin{1})`.
- Lines 58-59: Preserve location; implemented by `obj.isGPU=isa(varargin{1},'gpuArray')`.
- Lines 61-62: Get dimensions; implemented by `obj.numRows=int64(size(varargin{1},1))`.
- Lines 65-66: Get non-zeroes; implemented by `[row,col,val]=find(varargin{1})`.
- Lines 68-69: Make the object; implemented by `obj.row=int64(row(:))`.
- Lines 75-76: Complain and bomb out; implemented by `error('input cannot be converted into RCV sparse format.')`.
- Lines 82-83: Empty object with specified dimensions; implemented by `obj.numRows=int64(varargin{1})`.
- Lines 92-93: Build an object from scratch; implemented by `obj.row=int64(varargin{1}(:))`.
- Lines 99-102: Decide the location; implemented by `obj.isGPU=isa(varargin{1},'gpuArray')|| isa(varargin{2},'gpuArray')|| isa(varargin{3},'gpuArray')`.
- Lines 104-105: Upload to GPU; implemented by `if obj.isGPU`.
- Lines 113-114: Complain and bomb out; implemented by `error('incorrect number of input arguments.')`.

### Control flow inferred from the code

- Line 47: conditional branch on `nargin==1`.
- Line 50: conditional branch on `isa(varargin{1},'rcv')`.
- Line 105: conditional branch on `obj.isGPU`.

### Key state/data transformations

- Lines 53: computes `obj` using `obj=varargin{1}`.
- Lines 59: computes `obj.isGPU` using `obj.isGPU=isa(varargin{1},'gpuArray')`.
- Lines 62: computes `obj.numRows` using `obj.numRows=int64(size(varargin{1},1))`.
- Lines 63: computes `obj.numCols` using `obj.numCols=int64(size(varargin{1},2))`.
- Lines 66: computes `[row,col,val]` using `[row,col,val]=find(varargin{1})`.
- Lines 69: computes `obj.row` using `obj.row=int64(row(:))`.
- Lines 70: computes `obj.col` using `obj.col=int64(col(:))`.
- Lines 71: computes `obj.val` using `obj.val=double(val(:))`.

### Local helper functions

- Line 121: `isnumeric()` — `function answer=isnumeric(obj)`. Always yes
  - Representative operation: `answer=true()`.
- Line 129: `ismatrix()` — `function answer=ismatrix(obj)`. Always yes
  - Representative operation: `answer=true()`.
- Line 137: `isfloat()` — `function answer=isfloat(obj)`. Always yes
  - Representative operation: `answer=true()`.
- Line 149: `grumble()` — `function grumble(varargin)`.
  - Representative operation: `if nargin==1`.
  - Representative operation: `if (~isa(varargin{1},'rcv'))&&((~isnumeric(varargin{1}))||(~ismatrix(varargin{1})))`.

## Parameters / inputs

- M -a Matlab matrix
- dim1 -number of rows
- dim2 -number of columns
- R -row indices of non-zero entries
- C -column indices of non-zero entries
- V -values corresponding to entries in R and C

## Outputs

- obj -an RCV sparse matrix object

## Implementation structure

- Creates an RCV (row-column-value storage) sparse matrix. Syntax:
- obj=rcv(M)
- obj=rcv(dim1,dim2)
- obj=rcv(R,C,V,dim1,dim2)
- M -a Matlab matrix
- dim1 -number of rows
- dim2 -number of columns
- R -row indices of non-zero entries
- C -column indices of non-zero entries
- V -values corresponding to entries in R and C
- obj -an RCV sparse matrix object
- Check consistency

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `classdef()`, `int64()`, `double()`, `grumble()`, `ismatrix()`, `row()`, `col()`, `val()`, `gpuArray()`, `true()`, `isfloat()`, `isscalar()`, `any()`.
