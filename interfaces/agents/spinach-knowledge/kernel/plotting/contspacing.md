# kernel/plotting/contspacing.m

- Signature: `[all_conts,pos_conts,neg_conts]=...`

## Purpose

Non-linear adaptive contour spacing. Useful for NMR data where small cross-peaks must be adequately contoured next to large diagonal peaks.

## Physical / mathematical content

## Numerical / algorithmic content

## Syntax

```matlab
[all_conts,pos_conts,neg_conts]=...
contspacing(smax,smin,delta,k,signs,ncont)
```

## Parameters / inputs

- smax -global maximum intensity in the spectrum
- smin -global minimum intensity in the spectrum
- delta -minimum and maximum elevation (as a fraction of the
- total intensity) of the contours above the baseline.
- A good starting value is [0.02 0.2 0.02 0.2]. The
- first pair of numbers refers to the positive conto-
- urs and the second pair to the negative ones.
- k -a coefficient that controls the curvature of the contour
- spacing function: k=1 corresponds to linear spacing and
- k>1 bends the spacing curve to increase the sampling den-
- sity near the baseline. A reasonable value is 2.
- signs -can be set to 'positive', 'negative' or 'both' -this
- will cause the corresponding contours to be returned.
- ncont -the number of contours, a reasonable value is 20

## Outputs

- all_conts -all contour levels, a row vector
- pos_conts -positive contour levels, a row vector
- neg_conts -negative contour levels, a row vector
- Note: the following functions are used to get contour levels
- pos_conts=delta(2)*smax*linspace(0,1,ncont).^k+smax*delta(1);
- neg_conts=delta(2)*smin*linspace(0,1,ncont).^k+smin*delta(1);

## Implementation structure

- Non-linear adaptive contour spacing. Useful for NMR data where small
- cross-peaks must be adequately contoured next to large diagonal peaks.
- [all_conts,pos_conts,neg_conts]=...
- contspacing(smax,smin,delta,k,signs,ncont)
- smax -global maximum intensity in the spectrum
- smin -global minimum intensity in the spectrum
- delta -minimum and maximum elevation (as a fraction of the
- total intensity) of the contours above the baseline.
- A good starting value is [0.02 0.2 0.02 0.2]. The
- first pair of numbers refers to the positive conto-
- urs and the second pair to the negative ones.
- k -a coefficient that controls the curvature of the contour
