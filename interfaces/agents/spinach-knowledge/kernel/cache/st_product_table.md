# kernel/cache/st_product_table.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/kernel/cache/st_product_table.m`
- Signature: `[pt_left,pt_right]=st_product_table(nlevels)`
- Total lines: 94

## Purpose

Structure coefficient tables for single transition operators. Syntax: [pt_left,pt_right]=st_product_table(nlevels)

## Physical / mathematical content

- Cache-management utilities. These files maintain Spinach temporary or persistent cache state used to avoid repeated expensive construction of large operators or metadata.

## Numerical / algorithmic content

- The file contains an explicit `grumble(...)` validator, which is Spinach convention for front-loading dimension, type, and regime checks before expensive linear-algebra work begins.
- The file also defines local helper function(s): `grumble()`. This usually means the public entry point is supported by tightly coupled validation or helper logic kept private to the file.

## Code-derived implementation details

### Comment-guided execution stages

- Lines 32-33: Check consistency; implemented by `grumble(nlevels)`.
- Lines 35-36: Generate cache record name; implemented by `own_path=mfilename('fullpath')`.
- Lines 41-42: Check the cache; implemented by `if exist(table_file,'file')`.
- Lines 44-45: Lift data from the cache if the file is already available; implemented by `load(table_file,'pt_left','pt_right')`.
- Lines 49-50: Get ST operators; implemented by `B=sin_tran(nlevels)`.
- Lines 52-53: Preallocate the arrays; implemented by `pt_left= zeros(nlevels^2,nlevels^2,nlevels^2)`.
- Lines 56-57: Get the structure coefficients; implemented by `for k=1:nlevels^2`.
- Lines 61-62: Left product action; implemented by `pt_left(n,m,k)= hdot(B{k},B{n}*B{m})`.
- Lines 64-65: Right product action; implemented by `pt_right(n,m,k)=hdot(B{k},B{m}*B{n})`.

### Control flow inferred from the code

- Line 42: conditional branch on `exist(table_file,'file')`.
- Line 57: `for` loop over `k=1:nlevels^2`.
- Line 58: `for` loop over `m=1:nlevels^2`.
- Line 59: `for` loop over `n=1:nlevels^2`.

### Key state/data transformations

- Lines 36: computes `own_path` using `own_path=mfilename('fullpath')`.
- Lines 38-39: computes `table_file` using `table_file=[own_path 'st_product_table_' num2str(nlevels) '.mat']`.
- Lines 50: computes `B` using `B=sin_tran(nlevels)`.
- Lines 53: computes `pt_left` using `pt_left= zeros(nlevels^2,nlevels^2,nlevels^2)`.
- Lines 54: computes `pt_right` using `pt_right=zeros(nlevels^2,nlevels^2,nlevels^2)`.
- Lines 62: computes `pt_left(n,m,k)` using `pt_left(n,m,k)= hdot(B{k},B{n}*B{m})`.
- Lines 65: computes `pt_right(n,m,k)` using `pt_right(n,m,k)=hdot(B{k},B{m}*B{n})`.

### Local helper functions

- Line 82: `grumble()` — `function grumble(nlevels)`. Из толстых книг нельзя узнать ничего нового. Толстые книги -это кладбище, в котором погребены отслужившие
  - Representative operation: `if (~isnumeric(nlevels))||(~isreal(nlevels))|| (~isscalar(nlevels))||(nlevels<1)||(mod(nlevels,1)~=0)`.
  - Representative operation: `(~isscalar(nlevels))||(nlevels<1)||(mod(nlevels,1)~=0)`.

## Parameters / inputs

- nlevels -the number of energy levels in the system

## Outputs

- pt_left, pt_right -structure coefficients in
- the following conventions:
- S{n}*S{m}=...+pt_left(n,m,k)*S{k}+...
- S{m}*S{n}=...+pt_right(n,m,k)*S{k}+...
- corresponding to the expansion of the left and the right multiplica-
- tive action by S{n} on S{m} as given in Eq 7.18 of the first edition
- of IK's book (normalisation is missing in the book, that's a typo).
- Numbering translation between single and double index is given by
- kq2lin and lin2kq functions.
- Note: these are expensive tables, disk cache is used automatically.

## Implementation structure

- Structure coefficient tables for single transition operators. Syntax:
- [pt_left,pt_right]=st_product_table(nlevels)
- nlevels -the number of energy levels in the system
- pt_left, pt_right -structure coefficients in
- the following conventions:
- S{n}*S{m}=...+pt_left(n,m,k)*S{k}+...
- S{m}*S{n}=...+pt_right(n,m,k)*S{k}+...
- corresponding to the expansion of the left and the right multiplica-
- tive action by S{n} on S{m} as given in Eq 7.18 of the first edition
- of IK's book (normalisation is missing in the book, that's a typo).
- Numbering translation between single and double index is given by
- kq2lin and lin2kq functions.

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `grumble()`, `mfilename()`, `own_path()`, `num2str()`, `exist()`, `load()`, `sin_tran()`, `pt_left()`, `hdot()`, `pt_right()`, `save()`, `isscalar()`.
