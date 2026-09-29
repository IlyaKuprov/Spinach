# examples/nmr_liquids/inad_three_spin.m

- Signature: `inad_three_spin()`
- Source: [examples/nmr_liquids/inad_three_spin.m](https://github.com/IlyaKuprov/Spinach/blob/main/examples/nmr_liquids/inad_three_spin.m)

## Purpose

A one-dimensional liquid-state INADEQUATE example for a generic three-spin `13C` system. The model has one coupled pair and a third uncoupled carbon; the source estimates calculation time in seconds.

## Implementation

The field is `9.4` T. The three `13C` spin shifts are `10`, `15`, and `20` ppm, with `J(1,2)=55` Hz and no other nonzero pair coupling in the model. Spinach uses the full `sphten-liouv` basis with no approximation. The one-dimensional INADEQUATE simulation observes `13C`, sets the sequence coupling to `55` Hz, and has no decoupled channel. The source acquisition values are `offset=1800` and `sweep=5000` (units are not stated), with 4096 points and zero filling to 16384. It applies an exponential window with parameter 5, Fourier transforms and plots the real spectrum in ppm.
