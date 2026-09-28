# kernel/plotting/contspacing.m

- Signature: `[all_conts,pos_conts,neg_conts]=...`

## Purpose

Computes adaptive positive and/or negative contour levels for spectra, allowing small cross-peaks to be contoured alongside large diagonal peaks.

## Syntax

```matlab
[all_conts,pos_conts,neg_conts]=contspacing(smax,smin,delta,k,signs,ncont)
```

## Parameters / inputs

- `smax` — global maximum spectrum intensity.
- `smin` — global minimum spectrum intensity.
- `delta` — four contour fractions: `[positive_min positive_max negative_min negative_max]`; each lies in [0,1], with each minimum no greater than its maximum. A suggested value is `[0.02 0.2 0.02 0.2]`.
- `k` — positive integer curvature exponent; `k=1` gives linear spacing, while `k>1` increases sampling density near the baseline. A suggested value is 2.
- `signs` — `'positive'`, `'negative'`, or `'both'`.
- `ncont` — positive integer number of levels per requested sign; a suggested value is 20.

## Numerical / algorithmic content

For `t=linspace(0,1,ncont)`, positive levels are `smax*(delta(1)+(delta(2)-delta(1))*t.^k)`, and negative levels are `smin*(delta(3)+(delta(4)-delta(3))*t.^k)`. A sign-specific output is empty if that sign is not requested or its corresponding extremum does not have that sign. `all_conts` concatenates the negative levels in reverse order with the positive levels. Inputs are checked for finite real scalar extrema, valid fractions, a positive integer `k` and `ncont`, and one of the three supported `signs` values.

## Outputs

- `all_conts` — requested contour levels, negative levels first in ascending order, followed by positive levels.
- `pos_conts` — positive contour levels, or empty.
- `neg_conts` — negative contour levels, or empty.
