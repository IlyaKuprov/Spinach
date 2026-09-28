# kernel/cache/ist_product_table.m

- Signature: `[product_table_left,product_table_right]=ist_product_table(mult)`

## Purpose

Construct structure coefficient tables for the associative envelopes of su(mult) algebras, describing left and right multiplication in an irreducible spherical tensor (IST) basis.

## Physical / mathematical content

The coefficients follow the conventions:

- `T{n}*T{m}=...+product_table_left(n,m,k)*T{k}+...`
- `T{m}*T{n}=...+product_table_right(n,m,k)*T{k}+...`

These correspond to the IST expansion of the left and right multiplicative action by `T{n}` on `T{m}`, as given in Eq 7.18 of the first edition of IK's book. The normalisation is missing in the book; this is a typo. Numbering translation between single and double indices is provided by `lm2lin` and `lin2lm`.

## Numerical / algorithmic content

For multiplicity two (spin-1/2), the function returns hard-coded `4×4×4` left and right tables to avoid a disk-cache hit. For other multiplicities, it loads existing tables from a `.mat` cache record when available. Otherwise, it obtains IST operators from `irr_sph_ten(mult)`, computes their norms using `hdot`, and populates both `mult^2 × mult^2 × mult^2` tables from normalized operator products. Normalizing the operators separately in these products helps handle extreme norms. It then attempts to save the tables as a `-v7.3` cache record; failure to write produces a warning but does not prevent the tables from being returned.

## Parameters / inputs

- `mult` — multiplicity (number of energy levels) of the spin in question. It must be a real integer greater than 1.

## Outputs

- `product_table_left` — structure coefficients for left multiplication, under the convention above.
- `product_table_right` — structure coefficients for right multiplication, under the convention above.

## Implementation structure

The function validates `mult` with `grumble`, returns the hard-coded tables if `mult==2`, and otherwise checks for a cached record before calculating and optionally saving the tables. These tables are expensive to calculate, so disk caching is used automatically; the smallest and most frequently used case is hard-coded.

<https://spindynamics.org/wiki/index.php?title=ist_product_table.m>