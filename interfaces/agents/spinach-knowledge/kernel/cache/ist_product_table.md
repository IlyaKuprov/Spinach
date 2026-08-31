# kernel/cache/ist_product_table.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/kernel/cache/ist_product_table.m`
- Signature: `[product_table_left,product_table_right]=ist_product_table(mult)`
- Total lines: 156

## Purpose

Structure coefficient tables for the associative envelopes of su(mult) algebras. Syntax: [product_table_left,product_table_right]=ist_product_table(mult)

## Physical / mathematical content

- Cache-management utilities. These files maintain Spinach temporary or persistent cache state used to avoid repeated expensive construction of large operators or metadata.

## Numerical / algorithmic content

- The file contains an explicit `grumble(...)` validator, which is Spinach convention for front-loading dimension, type, and regime checks before expensive linear-algebra work begins.
- The file also defines local helper function(s): `grumble()`. This usually means the public entry point is supported by tightly coupled validation or helper logic kept private to the file.

## Code-derived implementation details

### Comment-guided execution stages

- Lines 37-38: Check consistency; implemented by `grumble(mult)`.
- Lines 40-42: Hard-coded for spin-1/2 to avoid expensive disk hit; implemented by `if mult==2`.
- Lines 44-45: Left side product table; implemented by `product_table_left=zeros(4,4,4)`.
- Lines 63-64: Right side product table; implemented by `product_table_right=zeros(4,4,4)`.
- Lines 84-85: Generate cache record name; implemented by `own_path=mfilename('fullpath')`.
- Lines 90-91: Check the cache; implemented by `if exist(table_file,'file')`.
- Lines 93-94: Lift data from the cache if the file is already available; implemented by `load(table_file,'product_table_left','product_table_right')`.
- Lines 98-99: Get IST operators; implemented by `T=irr_sph_ten(mult)`.
- Lines 101-102: Precompute norms; implemented by `norms=zeros(mult^2,1)`.
- Lines 107-108: Preallocate IST product tables; implemented by `product_table_left= zeros(mult^2,mult^2,mult^2)`.
- Lines 111-112: Populate product tables; implemented by `for k=1:mult^2`.
- Lines 116-119: Left product action: carefully tiptoe around extreme norms; implemented by `product_table_left(n,m,k)= norms(n)*hdot((T{k}/norms(k)), (T{n}/norms(n))* (T{m}/norms(m)))`.
- Lines 121-124: Right product action: carefully tiptoe around extreme norms; implemented by `product_table_right(n,m,k)=norms(n)*hdot((T{k}/norms(k)), (T{m}/norms(m))* (T{n}/norms(n)))`.

### Control flow inferred from the code

- Line 42: conditional branch on `mult==2`.
- Line 91: conditional branch on `exist(table_file,'file')`.
- Line 103: `for` loop over `n=1:mult^2`.
- Line 112: `for` loop over `k=1:mult^2`.
- Line 113: `for` loop over `m=1:mult^2`.
- Line 114: `for` loop over `n=1:mult^2`.

### Key state/data transformations

- Lines 45: computes `product_table_left` using `product_table_left=zeros(4,4,4)`.
- Lines 46: computes `product_table_left(:,:,1)` using `product_table_left(:,:,1)= [1.0 0 0 0`.
- Lines 50: computes `product_table_left(:,:,2)` using `product_table_left(:,:,2)= [ 0 1.0 0 0`.
- Lines 54: computes `product_table_left(:,:,3)` using `product_table_left(:,:,3)= [ 0 0 1.0 0`.
- Lines 58: computes `product_table_left(:,:,4)` using `product_table_left(:,:,4)= [ 0 0 0 1.0`.
- Lines 64: computes `product_table_right` using `product_table_right=zeros(4,4,4)`.
- Lines 65: computes `product_table_right(:,:,1)` using `product_table_right(:,:,1)=[1.0 0 0 0`.
- Lines 69: computes `product_table_right(:,:,2)` using `product_table_right(:,:,2)=[ 0 1.0 0 0`.
- Lines 73: computes `product_table_right(:,:,3)` using `product_table_right(:,:,3)=[ 0 0 1.0 0`.
- Lines 77: computes `product_table_right(:,:,4)` using `product_table_right(:,:,4)=[ 0 0 0 1.0`.
- Lines 85: computes `own_path` using `own_path=mfilename('fullpath')`.
- Lines 87-88: computes `table_file` using `table_file=[own_path 'ist_product_table_' num2str(mult) '.mat']`.
- Lines 99: computes `T` using `T=irr_sph_ten(mult)`.
- Lines 102: computes `norms` using `norms=zeros(mult^2,1)`.
- Lines 104: computes `norms(n)` using `norms(n)=sqrt(hdot(T{n},T{n}))`.
- Lines 117-119: computes `product_table_left(n,m,k)` using `product_table_left(n,m,k)= norms(n)*hdot((T{k}/norms(k)), (T{n}/norms(n))* (T{m}/norms(m)))`.
- Lines 122-124: computes `product_table_right(n,m,k)` using `product_table_right(n,m,k)=norms(n)*hdot((T{k}/norms(k)), (T{m}/norms(m))* (T{n}/norms(n)))`.

### Local helper functions

- Line 142: `grumble()` — `function grumble(mult)`. According to Oxford Chemistry folklore, Peter Atkins has once asked the following question at an interview for a Lecturer post: "What is it that
  - Representative operation: `if (numel(mult)~=1)||(~isnumeric(mult))|| (~isreal(mult))||(mult<2)||(mod(mult,1)~=0)`.
  - Representative operation: `(~isreal(mult))||(mult<2)||(mod(mult,1)~=0)`.

## Parameters / inputs

- mult -multiplicity (number of energy levels) of
- the spin in question

## Outputs

- product_table_left
- product_table_right -structure coefficients in
- the following conventions:
- T{n}*T{m}=...+product_table_left(n,m,k)*T{k}+...
- T{m}*T{n}=...+product_table_right(n,m,k)*T{k}+...
- corresponding to the IST expansion of the left and the right multiplica-
- tive action by T{n} on T{m} as given in Eq 7.18 of the first edition of
- IK's book (normalisation is missing in the book, that's a typo). Number-
- ing translation between single and double index is given by lm2lin and
- lin2lm functions.
- Note: these are expensive tables, disk cache is used automatically,
- the smallest and most frequently used case is hard-coded.

## Implementation structure

- Structure coefficient tables for the associative envelopes of su(mult)
- algebras. Syntax:
- [product_table_left,product_table_right]=ist_product_table(mult)
- mult -multiplicity (number of energy levels) of
- the spin in question
- product_table_left
- product_table_right -structure coefficients in
- the following conventions:
- T{n}*T{m}=...+product_table_left(n,m,k)*T{k}+...
- T{m}*T{n}=...+product_table_right(n,m,k)*T{k}+...
- corresponding to the IST expansion of the left and the right multiplica-
- tive action by T{n} on T{m} as given in Eq 7.18 of the first edition of

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `grumble()`, `product_table_left()`, `product_table_right()`, `mfilename()`, `own_path()`, `num2str()`, `exist()`, `load()`, `irr_sph_ten()`, `norms()`, `hdot()`, `save()`.
