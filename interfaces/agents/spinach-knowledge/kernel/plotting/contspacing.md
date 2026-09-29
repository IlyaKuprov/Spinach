# kernel/plotting/contspacing.m

- Source: [kernel/plotting/contspacing.m](https://github.com/IlyaKuprov/Spinach/blob/main/kernel/plotting/contspacing.m) · [Wiki](https://spindynamics.org/wiki/index.php?title=contspacing.m)
- Signature: `[all_conts,pos_conts,neg_conts]=contspacing(smax,smin,delta,k,signs,ncont)`

## Purpose

Builds nonlinear positive and/or negative contour levels so weak cross-peaks can be contoured alongside strong peaks. It returns levels only; it does not plot a spectrum.

## Inputs

- `smax`, `smin` — global maximum and minimum spectrum intensities.
- `delta` — four finite fractions in [0,1], ordered as [positive minimum, positive maximum, negative minimum, negative maximum]. Each pair must be ascending. The documented starting value is `[0.02 0.2 0.02 0.2]`.
- `k` — positive integer curvature exponent. `k=1` gives linear spacing; `k>1` makes the power curve nonlinear.
- `signs` — `'positive'`, `'negative'`, or `'both'`.
- `ncont` — positive integer number of levels on each requested side; 20 is a documented reasonable value.

## Level construction and ordering

For `u=linspace(0,1,ncont)`, positive levels run from `smax*delta(1)` to `smax*delta(2)` using `u.^k`. Negative levels run from `smin*delta(3)` to `smin*delta(4)`; because `smin` is negative, that second endpoint is more negative. With `k>1`, levels cluster toward the start of each progression (the less intense edge of that side). Positive levels are produced only when `smax>0`; negative levels only when `smin<0`; an absent or unrequested side is empty.

The returned arrays are row vectors. `pos_conts` is ascending from lower to higher positive intensity; `neg_conts` is reversed before concatenation, so `all_conts` places the most negative levels first, then the less negative levels, then the positive levels. Empty sides contribute no entries.

## Outputs

- `all_conts` — combined requested levels, ordered from negative to positive.
- `pos_conts` — positive levels, or an empty array.
- `neg_conts` — negative levels, or an empty array.
