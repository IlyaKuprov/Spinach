# kernel/utilities/dfpt.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/kernel/utilities/dfpt.m`
- Signature: `subgraphs=dfpt(conmatrix,max_sg_size)`
- Total lines: 110

## Purpose

Graph partitioning module. Analyzes the system connectivity graph and creates a list of all connected subgraphs of up to the user-specified size by crawling the graph in all available directions. Syntax: subgraphs=dfpt(conmatrix,max_sg_size)

## Physical / mathematical content

- General mathematical and infrastructure utilities. This area contains finite differences, perturbation theory, graph algorithms, spectral densities, tensor algebra, hash/report helpers, and other reusable numerical components.

## Numerical / algorithmic content

- The file contains an explicit `grumble(...)` validator, which is Spinach convention for front-loading dimension, type, and regime checks before expensive linear-algebra work begins.
- The file also defines local helper function(s): `size()`. This usually means the public entry point is supported by tightly coupled validation or helper logic kept private to the file.

## Code-derived implementation details

### Comment-guided execution stages

- Lines 28-29: Check consistency; implemented by `grumble(conmatrix,max_sg_size)`.
- Lines 31-32: Start at each spin; implemented by `nspins=size(conmatrix,2)`.
- Lines 35-36: Crawl the graph; implemented by `for sg_size=2:max_sg_size`.
- Lines 38-39: Preallocate grown set; implemented by `grown_set=cell(size(subgraphs,1),1)`.
- Lines 41-42: Loop over subgraphs; implemented by `for n=1:size(subgraphs,1)`.
- Lines 44-45: Get current subgraph; implemented by `subgraph=subgraphs(n,:)`.
- Lines 47-48: Find spins reachable from current subgraph; implemented by `neighbours=any(conmatrix(:,subgraph),2)`.
- Lines 52-53: Grow in every direction; implemented by `if isempty(neighbours)`.
- Lines 55-56: Isolated subgraphs get a dummy index; implemented by `grown_set{n}=[subgraph subgraph(end)]`.
- Lines 60-61: Subgraphs with neighbours get neighbours; implemented by `grown_set{n}=[repmat(subgraph,[numel(neighbours) 1]) neighbours]`.
- Lines 67-68: Merge grown set; implemented by `subgraphs=cell2mat(grown_set)`.
- Lines 74-75: Get sparse row index; implemented by `row_index=uint32(1:size(subgraphs,1))'`.
- Lines 79-80: Get sparse column index; implemented by `nsg=size(subgraphs,1); subgraphs=subgraphs(:)`.
- Lines 82-83: Return subgraph array as a sparse logical matrix; implemented by `subgraphs=sparse(row_index,subgraphs,1,nsg,nspins)`.

### Control flow inferred from the code

- Line 36: `for` loop over `sg_size=2:max_sg_size`.
- Line 42: `for` loop over `n=1:size(subgraphs,1)`.
- Line 53: conditional branch on `isempty(neighbours)`.

### Key state/data transformations

- Lines 32: computes `nspins` using `nspins=size(conmatrix,2)`.
- Lines 33: computes `subgraphs` using `subgraphs=uint32(1:nspins)'`.
- Lines 39: computes `grown_set` using `grown_set=cell(size(subgraphs,1),1)`.
- Lines 45: computes `subgraph` using `subgraph=subgraphs(n,:)`.
- Lines 48: computes `neighbours` using `neighbours=any(conmatrix(:,subgraph),2)`.
- Lines 49: computes `neighbours(subgraph)` using `neighbours(subgraph)=false()`.
- Lines 56: computes `grown_set{n}` using `grown_set{n}=[subgraph subgraph(end)]`.
- Lines 75: computes `row_index` using `row_index=uint32(1:size(subgraphs,1))'`.
- Lines 80: computes `nsg` using `nsg=size(subgraphs,1); subgraphs=subgraphs(:)`.

### Local helper functions

- Line 89: `grumble()` — `function grumble(conmatrix,max_sg_size)`.
  - Representative operation: `if ~islogical(conmatrix)`.
  - Representative operation: `error('conmatrix must be a logical matrix.')`.

## Parameters / inputs

- conmatrix -[nspins x nspins] matrix with 1 for
- connected spins and 0 elsewhere
- max_sg_size -maximum connected subgraph size

## Outputs

- subgraphs -[n_subgraphs x nspins] matrix; each
- row contains 1 for spins that belong
- to the subgraph and 0 for spins that
- do not

## Implementation structure

- Graph partitioning module. Analyzes the system connectivity graph and
- creates a list of all connected subgraphs of up to the user-specified
- size by crawling the graph in all available directions. Syntax:
- subgraphs=dfpt(conmatrix,max_sg_size)
- conmatrix -[nspins x nspins] matrix with 1 for
- connected spins and 0 elsewhere
- max_sg_size -maximum connected subgraph size
- subgraphs -[n_subgraphs x nspins] matrix; each
- row contains 1 for spins that belong
- to the subgraph and 0 for spins that
- do not
- Check consistency

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `grumble()`, `uint32()`, `subgraphs()`, `any()`, `conmatrix()`, `neighbours()`, `false()`, `subgraph()`, `cell2mat()`, `row_index()`, `logical()`, `islogical()`, `isscalar()`.
