# examples/nmr_liquids/cosy90_rotenone.m

- Signature: `cosy90_rotenone()`

## Purpose

COSY-90 spectrum of rotenone. Calculation time: minutes. Source assignment: [doi:10.1002/jhet.5570250160](http://dx.doi.org/10.1002/jhet.5570250160).

## Physical / mathematical content

This is a homonuclear liquid-state `1H` COSY-90 simulation for a 22-spin rotenone model. The source specifies proton chemical shifts and scalar couplings, then calls Spinach's liquid-state COSY simulation to generate a two-dimensional free-induction signal.

## Numerical / algorithmic content

The calculation sets field value 5.9 and uses the `sphten-liouv` formalism with IK-2 approximation, scalar-coupling connectivity, proximity level 1, the greedy option, and three S3 symmetry groups over spins 14–16, 17–19, and 20–22. It sets a 90-degree sequence angle, offset 1200, sweep 2000, 512 points and 2048 zero-fill points on both axes. Cosine apodisation precedes a shifted 2D FFT; the plotted spectrum is its real part.

## Implementation structure

The function explicitly defines 22 proton sites, their shifts and listed couplings, configures the basis and sequence, builds the Spinach system, and runs `liquid(...,@cosy,...,'nmr')`. It then apodises, transforms and plots the data.
