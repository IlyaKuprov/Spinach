# kernel/utilities/distrib_dim.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/kernel/utilities/distrib_dim.m`
- Signature: `A=distrib_dim(A,dim)`
- Total lines: 77

## Purpose

Distributes an array in the user-specified dimension for parallel processing using spmd. Syntax: A=distrib_dim(A,dim)

## Physical / mathematical content

- General mathematical and infrastructure utilities. This area contains finite differences, perturbation theory, graph algorithms, spectral densities, tensor algebra, hash/report helpers, and other reusable numerical components.

## Numerical / algorithmic content

- The implementation explicitly addresses performance engineering through parallel or GPU execution, which matters because Spinach operators can become extremely large after basis expansion or powder/spatial lifting.
- The file contains an explicit `grumble(...)` validator, which is Spinach convention for front-loading dimension, type, and regime checks before expensive linear-algebra work begins.
- The file also defines local helper function(s): `grumble()`. This usually means the public entry point is supported by tightly coupled validation or helper logic kept private to the file.

## Code-derived implementation details

### Comment-guided execution stages

- Lines 23-24: Check consistency; implemented by `grumble(A,dim)`.
- Lines 26-27: Get the size; implemented by `size_A=size(A)`.
- Lines 29-30: Set the stage; implemented by `spmd`.
- Lines 32-33: Codistributor with default partitioning; implemented by `defpart=codistributor1d.unsetPartition`.
- Lines 36-37: Start local parts; implemented by `LocalParts=1`.
- Lines 41-42: Retrieve codistribution index ranges; implemented by `CoD=CoD{1}; partLimits=[0 cumsum(CoD.Partition)]`.
- Lines 44-45: Dimension-agnostic index array; implemented by `idx=repmat({':'},1,ndims(A))`.
- Lines 47-48: Pull local parts; implemented by `for n=1:numel(LocalParts)`.
- Lines 53-54: Build a distributed array; implemented by `A=distributed(LocalParts,dim)`.

### Control flow inferred from the code

- Line 48: `for` loop over `n=1:numel(LocalParts)`.

### Key state/data transformations

- Lines 27: computes `size_A` using `size_A=size(A)`.
- Lines 33: computes `defpart` using `defpart=codistributor1d.unsetPartition`.
- Lines 34: computes `CoD` using `CoD=codistributor1d(dim,defpart,size_A)`.
- Lines 37: computes `LocalParts` using `LocalParts=1`.
- Lines 45: computes `idx` using `idx=repmat({':'},1,ndims(A))`.
- Lines 49: computes `idx{dim}` using `idx{dim}=partLimits(n)+1:partLimits(n+1)`.
- Lines 50: computes `LocalParts{n}` using `LocalParts{n}=A(idx{:})`.
- Lines 54: computes `A` using `A=distributed(LocalParts,dim)`.

### Local helper functions

- Line 59: `grumble()` — `function grumble(A,dim)`. The hostages decided amongst themselves that the two to
  - Representative operation: `if ~isnumeric(A), error('A must be numeric.'); end`.
  - Representative operation: `if (~isreal(dim))||(~isscalar(dim))|| (mod(dim,1)~=0)||(dim<1)`.

## Parameters / inputs

- A -a numerical array
- dim -the distribution dimension
- Output:
- A -a distributed numerical array
- Mathworks, Inc.

## Implementation structure

- Distributes an array in the user-specified dimension
- for parallel processing using spmd. Syntax:
- A=distrib_dim(A,dim)
- A -a numerical array
- dim -the distribution dimension
- Output:
- A -a distributed numerical array
- Mathworks, Inc.
- Check consistency
- Get the size
- Set the stage
- Codistributor with default partitioning

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `grumble()`, `codistributor1d()`, `cumsum()`, `ndims()`, `partLimits()`, `distributed()`, `isscalar()`.
