# kernel/plotting/contspacing.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/kernel/plotting/contspacing.m`
- Signature: `[all_conts,pos_conts,neg_conts]=...`
- Total lines: 108

## Purpose

Non-linear adaptive contour spacing. Useful for NMR data where small cross-peaks must be adequately contoured next to large diagonal peaks.

## Physical / mathematical content

- This file belongs to the `kernel` part of Spinach. Its role should be read together with nearby files in the same directory, which usually share a common physical regime or infrastructure purpose.

## Numerical / algorithmic content

- The file contains an explicit `grumble(...)` validator, which is Spinach convention for front-loading dimension, type, and regime checks before expensive linear-algebra work begins.

## Code-derived implementation details

### Comment-guided execution stages

- Lines 50-51: Check consistency; implemented by `grumble(smax,smin,delta,k,signs,ncont)`.
- Lines 53-54: Compute positive contour levels; implemented by `if (smax>0)&&(strcmp(signs,'positive')||strcmp(signs,'both'))`.
- Lines 60-61: Compute negative contour levels; implemented by `if (smin<0)&&(strcmp(signs,'negative')||strcmp(signs,'both'))`.
- Lines 67-68: Merge contour level arrays; implemented by `all_conts=[neg_conts(end:-1:1) pos_conts]`.

### Control flow inferred from the code

- Line 54: conditional branch on `(smax>0)&&(strcmp(signs,'positive')||strcmp(signs,'both'))`.
- Line 61: conditional branch on `(smin<0)&&(strcmp(signs,'negative')||strcmp(signs,'both'))`.

### Key state/data transformations

- Lines 55: computes `pos_conts` using `pos_conts=(delta(2)-delta(1))*smax*linspace(0,1,ncont).^k+smax*delta(1)`.
- Lines 62: computes `neg_conts` using `neg_conts=(delta(4)-delta(3))*smin*linspace(0,1,ncont).^k+smin*delta(3)`.
- Lines 68: computes `all_conts` using `all_conts=[neg_conts(end:-1:1) pos_conts]`.

### Local helper functions

- Line 73: `grumble()` — `function grumble(smax,smin,delta,k,signs,ncont)`.
  - Representative operation: `if (~isnumeric(smax))||(~isscalar(smax))||(~isreal(smax))||(~isfinite(smax))`.
  - Representative operation: `error('smax must be a finite real scalar.')`.

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

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `contspacing()`, `grumble()`, `strcmp()`, `delta()`, `neg_conts()`, `isscalar()`, `any()`, `ischar()`, `ismember()`.
