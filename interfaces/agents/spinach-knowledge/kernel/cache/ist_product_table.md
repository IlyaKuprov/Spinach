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
