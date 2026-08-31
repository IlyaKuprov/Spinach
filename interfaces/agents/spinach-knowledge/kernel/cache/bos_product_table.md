# kernel/cache/bos_product_table.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/kernel/cache/bos_product_table.m`
- Signature: `[product_table_left,...`
- Total lines: 111

## Purpose

Structure coefficient tables for the associative envelopes of truncated Weyl algebras spanned by orthogonalised bosonic mo- nomials. Syntax: [product_table_left,... product_table_right]=bos_product_table(nlevels)

## Physical / mathematical content

- Cache-management utilities. These files maintain Spinach temporary or persistent cache state used to avoid repeated expensive construction of large operators or metadata.

## Numerical / algorithmic content

- The file contains an explicit `grumble(...)` validator, which is Spinach convention for front-loading dimension, type, and regime checks before expensive linear-algebra work begins.
- The file also defines local helper function(s): `grumble()`. This usually means the public entry point is supported by tightly coupled validation or helper logic kept private to the file.

## Code-derived implementation details

### Comment-guided execution stages

- Lines 38-39: Check consistency; implemented by `grumble(nlevels)`.
- Lines 41-42: Generate cache record name; implemented by `own_path=mfilename('fullpath')`.
- Lines 47-48: Check the cache; implemented by `if exist(table_file,'file')`.
- Lines 50-51: Lift data from the cache if the file is already available; implemented by `load(table_file,'product_table_left','product_table_right')`.
- Lines 55-56: Orthogonalised bosonic monomials; implemented by `B=boson_ortho(nlevels)`.
- Lines 58-59: Precompute norms; implemented by `norms=zeros(nlevels^2,1)`.
- Lines 64-65: Preallocate OBM product tables; implemented by `product_table_left= zeros(nlevels^2,nlevels^2,nlevels^2)`.
- Lines 68-69: Populate product tables; implemented by `for k=1:nlevels^2`.
- Lines 73-76: Left product action; implemented by `product_table_left(n,m,k)= norms(n)*hdot((B{k}/norms(k)), (B{n}/norms(n))* (B{m}/norms(m)))`.
- Lines 78-81: Right product action; implemented by `product_table_right(n,m,k)=norms(n)*hdot((B{k}/norms(k)), (B{m}/norms(m))* (B{n}/norms(n)))`.

### Control flow inferred from the code

- Line 48: conditional branch on `exist(table_file,'file')`.
- Line 60: `for` loop over `n=1:nlevels^2`.
- Line 69: `for` loop over `k=1:nlevels^2`.
- Line 70: `for` loop over `m=1:nlevels^2`.
- Line 71: `for` loop over `n=1:nlevels^2`.

### Key state/data transformations

- Lines 36: computes `product_table_right]` using `product_table_right]=bos_product_table(nlevels)`.
- Lines 42: computes `own_path` using `own_path=mfilename('fullpath')`.
- Lines 44-45: computes `table_file` using `table_file=[own_path 'bos_product_table_' num2str(nlevels) '.mat']`.
- Lines 56: computes `B` using `B=boson_ortho(nlevels)`.
- Lines 59: computes `norms` using `norms=zeros(nlevels^2,1)`.
- Lines 61: computes `norms(n)` using `norms(n)=sqrt(hdot(B{n},B{n}))`.
- Lines 65: computes `product_table_left` using `product_table_left= zeros(nlevels^2,nlevels^2,nlevels^2)`.
- Lines 66: computes `product_table_right` using `product_table_right=zeros(nlevels^2,nlevels^2,nlevels^2)`.
- Lines 74-76: computes `product_table_left(n,m,k)` using `product_table_left(n,m,k)= norms(n)*hdot((B{k}/norms(k)), (B{n}/norms(n))* (B{m}/norms(m)))`.
- Lines 79-81: computes `product_table_right(n,m,k)` using `product_table_right(n,m,k)=norms(n)*hdot((B{k}/norms(k)), (B{m}/norms(m))* (B{n}/norms(n)))`.

### Local helper functions

- Line 99: `grumble()` — `function grumble(nlevels)`. Rage Against the Machine [a rock band] never specified what type of machine they were furious with, but I rec-
  - Representative operation: `if (~isnumeric(nlevels))||(~isreal(nlevels))|| (~isscalar(nlevels))||(nlevels<1)||(mod(nlevels,1)~=0)`.
  - Representative operation: `(~isscalar(nlevels))||(nlevels<1)||(mod(nlevels,1)~=0)`.

## Parameters / inputs

- nlevels -number of bosonic ladder population levels

## Outputs

- product_table_left
- product_table_right -structure coefficients in
- the following conventions:
- B{n}*B{m}=...+product_table_left(n,m,k)*B{k}+...
- B{m}*B{n}=...+product_table_right(n,m,k)*B{k}+...
- corresponding to the expansion of the left and the right mul-
- tiplicative action by B{n} on B{m} as given in Eq 7.18 of the
- first edition of IK's book (normalisation is missing in the
- book, that's a typo).
- Note: these are expensive tables, a disk cache is created and
- used automatically.

## Implementation structure

- Structure coefficient tables for the associative envelopes of
- truncated Weyl algebras spanned by orthogonalised bosonic mo-
- nomials. Syntax:
- [product_table_left,...
- product_table_right]=bos_product_table(nlevels)
- nlevels -number of bosonic ladder population levels
- product_table_left
- product_table_right -structure coefficients in
- the following conventions:
- B{n}*B{m}=...+product_table_left(n,m,k)*B{k}+...
- B{m}*B{n}=...+product_table_right(n,m,k)*B{k}+...
- corresponding to the expansion of the left and the right mul-

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `bos_product_table()`, `grumble()`, `mfilename()`, `own_path()`, `num2str()`, `exist()`, `load()`, `boson_ortho()`, `norms()`, `hdot()`, `product_table_left()`, `product_table_right()`, `save()`, `isscalar()`.
