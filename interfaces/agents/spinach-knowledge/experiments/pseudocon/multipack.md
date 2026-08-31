# experiments/pseudocon/multipack.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/experiments/pseudocon/multipack.m`
- Signature: `Ilm=multipack(ranks,moments)`
- Total lines: 66

## Purpose

Packs multipole moments from a linear stream into a cell array that is arranged by ranks. Syntax: Ilm=multipack(ranks,moments)

## Physical / mathematical content

- Paramagnetic-pseudocontact inference routines. The mathematics includes inverse problems, tensor parameterisation, interpolation, and regularisation.

## Numerical / algorithmic content

- The file contains an explicit `grumble(...)` validator, which is Spinach convention for front-loading dimension, type, and regime checks before expensive linear-algebra work begins.
- The file also defines local helper function(s): `grumble()`. This usually means the public entry point is supported by tightly coupled validation or helper logic kept private to the file.

## Code-derived implementation details

### Comment-guided execution stages

- Lines 33-34: Check consistency; implemented by `grumble(ranks,moments)`.
- Lines 36-37: Set the cell array dimensions; implemented by `Ilm=cell(size(ranks))`.
- Lines 39-40: Unpack the ranks; implemented by `current_position=0`.

### Control flow inferred from the code

- Line 41: `for` loop over `k=1:numel(ranks)`.

### Key state/data transformations

- Lines 37: computes `Ilm` using `Ilm=cell(size(ranks))`.
- Lines 40: computes `current_position` using `current_position=0`.
- Lines 42: computes `Ilm{k}` using `Ilm{k}=moments((current_position+1):(current_position+2*ranks(k)+1))`.

### Local helper functions

- Line 49: `grumble()` — `function grumble(ranks,moments)`.
  - Representative operation: `if (~isnumeric(ranks))||(~isreal(ranks))|| (~isrow(ranks))||any(mod(ranks,1)~=0)|| (numel(unique(ranks))~=numel(ranks))||any(ranks<0)`.
  - Representative operation: `(~isrow(ranks))||any(mod(ranks,1)~=0)|| (numel(unique(ranks))~=numel(ranks))||any(ranks<0)`.

## Parameters / inputs

- ranks -a vector of spherical ranks present, e.g. [0 1 2]
- moments -a vector of multipole moments for each rank,
- arranged in a linear stream

## Outputs

- Ilm -a cell array of vector corresponding to the multipole
- moments defined in http://dx.doi.org/10.1039/c6cp05437d
- for L=0, one element
- for L=1, three elements
- for L=2, five elements
- et cetera.

## Implementation structure

- Packs multipole moments from a linear stream into a cell array
- that is arranged by ranks. Syntax:
- Ilm=multipack(ranks,moments)
- ranks -a vector of spherical ranks present, e.g. [0 1 2]
- moments -a vector of multipole moments for each rank,
- arranged in a linear stream
- Ilm -a cell array of vector corresponding to the multipole
- moments defined in http://dx.doi.org/10.1039/c6cp05437d
- for L=0, one element
- for L=1, three elements
- for L=2, five elements
- et cetera.

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `grumble()`, `moments()`, `ranks()`, `isrow()`, `any()`.
