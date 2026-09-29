# examples/nmr_liquids/pa_sucrose.m

- Signature: `pa_sucrose()`
- Source: [examples/nmr_liquids/pa_sucrose.m](https://github.com/IlyaKuprov/Spinach/blob/main/examples/nmr_liquids/pa_sucrose.m)

## What it models and loads

A simulated liquid-state 1H pulse-acquire spectrum of sucrose. The wrapper parses `../standard_systems/sucrose.log`, described in its source comment as a vacuum DFT calculation, and passes the parsed data to `g2spinach` for 1H spin-system properties. This is a computational-chemistry input, not an experimental spectrum. The source comment estimates a run time of seconds; that estimate was not measured here.

## Spin system and relaxation

The wrapper sets `options.min_j=1.0`, passes 31.8 as an additional `g2spinach` argument (its meaning and units are not named in this wrapper), and sets the magnetic field to 14.1 (unit not stated). The basis is spherical-tensor Liouville space, IK-2 approximation, scalar-coupling connectivity and proximity level 2. Redfield relaxation is enabled, equilibrium is zero, retained terms are `secular`, and the correlation time is 1e-9 s (1 ns). The proximity cutoff is 4.0; its unit is not stated.

## Acquisition and processing

The initial density operator and receiver are the 1H raising state, with no decoupling. The wrapper calls `liquid(...,@acquire,...,'nmr')`; acquisition uses offset 1800, sweep 5000 and 8,192 points, with zero filling to 65,536. Offset and sweep are not unit-labelled in the wrapper. The displayed axis is ppm and inverted. It applies exponential apodisation with parameter 6, Fourier transforms the FID and plots the real spectrum.

## Output and limits

The source produces a plotted simulation, not a measured spectrum or saved data file. It does not report numerical peaks or relaxation rates. The positional `g2spinach` argument 31.8 is recorded without assigning it a meaning the wrapper does not specify; details of the parser and acquisition callback are likewise outside this file.
