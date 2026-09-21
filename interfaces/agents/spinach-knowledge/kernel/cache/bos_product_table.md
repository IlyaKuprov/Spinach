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
