# kernel/utilities/polinfo.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/kernel/utilities/polinfo.m`
- Signature: `polinfo(p,level,label)`
- Total lines: 111

## Purpose

Draws an ASCII diagram of a polyadic object. Syntax: polinfo(p)

## Physical / mathematical content

- General mathematical and infrastructure utilities. This area contains finite differences, perturbation theory, graph algorithms, spectral densities, tensor algebra, hash/report helpers, and other reusable numerical components.

## Numerical / algorithmic content

- The file contains an explicit `grumble(...)` validator, which is Spinach convention for front-loading dimension, type, and regime checks before expensive linear-algebra work begins.
- The file also defines local helper function(s): `grumble()`. This usually means the public entry point is supported by tightly coupled validation or helper logic kept private to the file.

## Code-derived implementation details

### Comment-guided execution stages

- Lines 22-23: Default settings; implemented by `if ~exist('level','var'), level=0; end`.
- Lines 26-27: Check consistency; implemented by `grumble(p,level,label)`.
- Lines 29-30: Set indentation level; implemented by `indent=repmat(' ',1,4*level)`.
- Lines 32-33: Get the dimensions; implemented by `[nrows,ncols]=size(p)`.
- Lines 35-36: Print label and size; implemented by `fprintf('%s%s [%dx%d]\n',indent,label,nrows,ncols)`.
- Lines 38-39: Print prefixes; implemented by `if isa(p,'polyadic')`.
- Lines 56-57: Print kron terms; implemented by `for n=1:numel(p.cores)`.
- Lines 72-73: Print suffixes; implemented by `if isa(p,'polyadic')`.

### Control flow inferred from the code

- Line 23: conditional branch on `~exist('level','var'), level=0; end`.
- Line 24: conditional branch on `~exist('label','var'), label='polyadic'; end`.
- Line 39: conditional branch on `isa(p,'polyadic')`.
- Line 40: conditional branch on `~isempty(p.prefix)`.
- Line 43: `for` loop over `n=1:numel(p.prefix)`.
- Line 44: conditional branch on `isa(p.prefix{n},'polyadic')`.
- Line 57: `for` loop over `n=1:numel(p.cores)`.
- Line 59: `for` loop over `k=1:numel(p.cores{n})`.
- Line 60: conditional branch on `isa(p.cores{n}{k},'polyadic')`.
- Line 73: conditional branch on `isa(p,'polyadic')`.
- Line 74: conditional branch on `~isempty(p.suffix)`.
- Line 77: `for` loop over `n=1:numel(p.suffix)`.
- Line 78: conditional branch on `isa(p.suffix{n},'polyadic')`.

### Key state/data transformations

- Lines 30: computes `indent` using `indent=repmat(' ',1,4*level)`.
- Lines 33: computes `[nrows,ncols]` using `[nrows,ncols]=size(p)`.

### Local helper functions

- Line 93: `grumble()` — `function grumble(p,level,label)`.
  - Representative operation: `if ~isa(p,'polyadic')`.
  - Representative operation: `error('p must be a polyadic object.')`.

## Parameters / inputs

- p -polyadic object

## Outputs

- an ASCII diagram to the console
- Note: polyadic objects can be huge, the code below
- avoids making memory copies.

## Implementation structure

- Draws an ASCII diagram of a polyadic object. Syntax:
- polinfo(p)
- p - polyadic object
- an ASCII diagram to the console
- Note: polyadic objects can be huge, the code below
- avoids making memory copies.
- Default settings
- Check consistency
- Set indentation level
- Get the dimensions
- Print label and size
- Print prefixes

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `exist()`, `grumble()`, `isscalar()`, `ischar()`.
