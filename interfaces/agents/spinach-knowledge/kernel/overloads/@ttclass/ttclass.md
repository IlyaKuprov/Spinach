# kernel/overloads/@ttclass/ttclass.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/kernel/overloads/@ttclass/ttclass.m`
- Signature: `tt=ttclass(coeff,kronterms,tolerance)`
- Total lines: 125

## Purpose

Creates an object of a tensor train class. A tensor train is a type of un-opened Kronecker product that behaves as a matrix or a vector of a very large dimension, but takes a reasonable amount of memory to store. See https://doi.org/10.1137/090752286 for further informa- tion. Syntax: tt=ttclass(coeff,kronterms,tolerance)

## Physical / mathematical content

- Tensor-train linear algebra. These files implement compressed high-dimensional operators and AMEn/SVD-based algebra in tensor-train format.

## Numerical / algorithmic content

- The file contains an explicit `grumble(...)` validator, which is Spinach convention for front-loading dimension, type, and regime checks before expensive linear-algebra work begins.
- The file also defines local helper function(s): `ncores()`, `ntrains()`, `numArgumentsFromSubscript()`, `grumble()`. This usually means the public entry point is supported by tightly coupled validation or helper logic kept private to the file.

## Code-derived implementation details

### Comment-guided execution stages

- Lines 57-58: Check the number of inputs; implemented by `if nargin==0, return`.
- Lines 63-64: Check consistency; implemented by `grumble(coeff,kronterms,tolerance)`.
- Lines 66-67: Convert Kronecker products into tensor trains; implemented by `tt.cores=kronterms`.
- Lines 74-75: Store the coefficients; implemented by `tt.coeff=coeff`.
- Lines 77-78: Store the tolerances; implemented by `tt.tolerance=tolerance`.
- Lines 80-81: Do not print diagnostics; implemented by `tt.debuglevel=0`.

### Control flow inferred from the code

- Line 58: conditional branch on `nargin==0, return`.
- Line 68: `for` loop over `n=1:size(kronterms,2)`.
- Line 69: `for` loop over `d=1:size(kronterms,1)`.

### Key state/data transformations

- Lines 67: computes `tt.cores` using `tt.cores=kronterms`.
- Lines 70: computes `tt.cores{d,n}` using `tt.cores{d,n}=reshape(full(kronterms{d,n}),[1 size(kronterms{d,n}) 1])`.
- Lines 75: computes `tt.coeff` using `tt.coeff=coeff`.
- Lines 78: computes `tt.tolerance` using `tt.tolerance=tolerance`.
- Lines 81: computes `tt.debuglevel` using `tt.debuglevel=0`.

### Local helper functions

- Line 86: `ncores()` — `function answer=ncores(obj)`. Number of trains in a sum
  - Representative operation: `answer=size(obj.cores,1)`.
- Line 91: `ntrains()` — `function answer=ntrains(obj)`. Matlab subsref bug workaround
  - Representative operation: `answer=size(obj.cores,2)`.
- Line 96: `numArgumentsFromSubscript()` — `function n=numArgumentsFromSubscript(obj,s,ic)`. Consistency enforcement
  - Representative operation: `n=builtin('numArgumentsFromSubscript',obj,s,ic)`.
- Line 105: `grumble()` — `function grumble(coeff,cores,tolerance)`.
  - Representative operation: `if (~isnumeric(coeff))||(~isrow(coeff))`.
  - Representative operation: `error('coeff must be a complex row vector.')`.

## Parameters / inputs

- coeff -coefficient in front of the spin operator,
- usually the interaction magnitude
- kronterms -column cell array of matrices whose Krone-
- cker product makes up the spin operator
- tolerance -maximum deviation in the 2-norm between the
- TT representation and the flat matrix repre-
- sentation that the TT format is allowed to
- introduce

## Outputs

- tt -tensor train object

## Header notes

- 1. If multiple columns are supplied in kronterms, multiple coeffi-
- cients are given in coeff, and multiple tolerances are given in
- tolerance, the resulting tensor train is assumed to be the sum
- of the individual tensor trains specified in different columns.
- 2. Tensor trains are exotic and capricious structures, do not use
- them unless you know what you are doing.

## Implementation structure

- Creates an object of a tensor train class. A tensor train is a type
- of un-opened Kronecker product that behaves as a matrix or a vector
- of a very large dimension, but takes a reasonable amount of memory
- to store. See https://doi.org/10.1137/090752286 for further informa-
- tion. Syntax:
- tt=ttclass(coeff,kronterms,tolerance)
- coeff -coefficient in front of the spin operator,
- usually the interaction magnitude
- kronterms -column cell array of matrices whose Krone-
- cker product makes up the spin operator
- tolerance -maximum deviation in the 2-norm between the
- TT representation and the flat matrix repre-

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `elseif()`, `grumble()`, `ncores()`, `ntrains()`, `numArgumentsFromSubscript()`, `builtin()`, `isrow()`, `iscell()`, `any()`, `cellfun()`, `cores()`.
