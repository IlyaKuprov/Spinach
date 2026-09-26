# kernel/cache/st_product_table.m

- Signature: `[pt_left,pt_right]=st_product_table(nlevels)`

## Purpose

Structure coefficient tables for single transition operators. Syntax: [pt_left,pt_right]=st_product_table(nlevels)

## Physical / mathematical content

- Cache-management utilities. These files maintain Spinach temporary or persistent cache state used to avoid repeated expensive construction of large operators or metadata.

## Numerical / algorithmic content

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
