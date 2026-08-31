# kernel/utilities/isoswap.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/kernel/utilities/isoswap.m`
- Signature: `[sys,inter]=isoswap(sys,inter,spins,new_iso)`
- Total lines: 116

## Purpose

Makes isotope replacements in the input structures. All interactions are automatically scaled as necessary. Syntax: [sys,inter]=isoswap(sys,inter,spins,new_iso)

## Physical / mathematical content

- General mathematical and infrastructure utilities. This area contains finite differences, perturbation theory, graph algorithms, spectral densities, tensor algebra, hash/report helpers, and other reusable numerical components.

## Numerical / algorithmic content

- The file contains an explicit `grumble(...)` validator, which is Spinach convention for front-loading dimension, type, and regime checks before expensive linear-algebra work begins.
- The file also defines local helper function(s): `grumble()`. This usually means the public entry point is supported by tightly coupled validation or helper logic kept private to the file.

## Code-derived implementation details

### Comment-guided execution stages

- Lines 29-30: Grumbler missing; implemented by `grumble(sys,inter,spins,new_iso)`.
- Lines 32-33: Wipe quadratic couplings -eigensystem representation; implemented by `if isfield(inter,'coupling')&&isfield(inter.coupling,'eigs')`.
- Lines 43-44: Wipe quadratic couplings -matrix representation; implemented by `if isfield(inter,'coupling')&&isfield(inter.coupling,'matrix')`.
- Lines 54-55: Wipe high-rank couplings; implemented by `if isfield(inter,'giant')`.
- Lines 65-66: Loop over spins; implemented by `for n=spins(:)'`.
- Lines 68-69: Get the gamma ratio; implemented by `gamma_ratio=spin(new_iso)/spin(sys.isotopes{n})`.
- Lines 71-72: Loop over coupling partners; implemented by `for k=setdiff(1:numel(sys.isotopes),n)`.
- Lines 74-75: Scale eigensystem representation; implemented by `if isfield(inter.coupling,'eigs')`.
- Lines 80-81: Scale matrix representation; implemented by `if isfield(inter.coupling,'matrix')`.
- Lines 88-89: Replace the isotope string; implemented by `sys.isotopes{n}=new_iso`.

### Control flow inferred from the code

- Line 33: conditional branch on `isfield(inter,'coupling')&&isfield(inter.coupling,'eigs')`.
- Line 34: `for` loop over `n=spins(:)'`.
- Line 35: conditional branch on `~isempty(inter.coupling.eigs{n,n})`.
- Line 44: conditional branch on `isfield(inter,'coupling')&&isfield(inter.coupling,'matrix')`.
- Line 45: `for` loop over `n=spins(:)'`.
- Line 46: conditional branch on `~isempty(inter.coupling.matrix{n,n})`.
- Line 55: conditional branch on `isfield(inter,'giant')`.
- Line 56: `for` loop over `n=spins(:)'`.
- Line 57: conditional branch on `~isempty(inter.giant.coef{n})`.
- Line 66: `for` loop over `n=spins(:)'`.
- Line 72: `for` loop over `k=setdiff(1:numel(sys.isotopes),n)`.
- Line 75: conditional branch on `isfield(inter.coupling,'eigs')`.
- Line 81: conditional branch on `isfield(inter.coupling,'matrix')`.

### Key state/data transformations

- Lines 36: computes `inter.coupling.eigs{n,n}` using `inter.coupling.eigs{n,n}=[]; inter.coupling.euler{n,n}=[]`.
- Lines 47: computes `inter.coupling.matrix{n,n}` using `inter.coupling.matrix{n,n}=[]`.
- Lines 58: computes `inter.giant.coef{n}` using `inter.giant.coef{n}={}; inter.giant.euler{n}={}`.
- Lines 69: computes `gamma_ratio` using `gamma_ratio=spin(new_iso)/spin(sys.isotopes{n})`.
- Lines 76: computes `inter.coupling.eigs{n,k}` using `inter.coupling.eigs{n,k}=gamma_ratio*inter.coupling.eigs{n,k}`.
- Lines 77: computes `inter.coupling.eigs{k,n}` using `inter.coupling.eigs{k,n}=gamma_ratio*inter.coupling.eigs{k,n}`.
- Lines 82: computes `inter.coupling.matrix{n,k}` using `inter.coupling.matrix{n,k}=gamma_ratio*inter.coupling.matrix{n,k}`.
- Lines 83: computes `inter.coupling.matrix{k,n}` using `inter.coupling.matrix{k,n}=gamma_ratio*inter.coupling.matrix{k,n}`.
- Lines 89: computes `sys.isotopes{n}` using `sys.isotopes{n}=new_iso`.

### Local helper functions

- Line 96: `grumble()` — `function grumble(sys,inter,spins,new_iso)`.
  - Representative operation: `if ~isfield(sys,'isotopes')`.
  - Representative operation: `error('isotope information missing from sys structure.')`.

## Parameters / inputs

- sys, inter -Spinach input data structures
- spins -a vector of integers specifying
- spin numbers
- new_iso -character string specifying
- the new isotope
- Output:
- sys, inter -Spinach input data structures
- Note: quadratic and higher order couplings are wiped and a warning
- is printed -those are not transferable.

## Implementation structure

- Makes isotope replacements in the input structures. All interactions
- are automatically scaled as necessary. Syntax:
- [sys,inter]=isoswap(sys,inter,spins,new_iso)
- sys, inter -Spinach input data structures
- spins -a vector of integers specifying
- spin numbers
- new_iso -character string specifying
- the new isotope
- Output:
- Note: quadratic and higher order couplings are wiped and a warning
- is printed -those are not transferable.
- Grumbler missing

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `grumble()`, `isfield()`, `spins()`, `int2str()`, `wiped()`, `spin()`, `setdiff()`, `isstruct()`, `isvector()`, `any()`, `ischar()`.
