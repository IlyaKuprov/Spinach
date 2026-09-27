# examples/nmr_liquids/ct_cosy_rotenone.m

- Signature: `ct_cosy_rotenone()`

## Purpose

CT COSY spectrum of rotenone using the assignment reported in [doi:10.1002/jhet.5570250160](http://dx.doi.org/10.1002/jhet.5570250160). Calculation time: minutes

## Physical / mathematical content

This is a constant-time homonuclear COSY simulation using the 22-proton rotenone assignment. The source specifies proton shifts and scalar couplings, then generates a two-dimensional signal with Spinach's liquid-state `ct_cosy` sequence.

## Numerical / algorithmic content

The model sets field value 5.9 and uses the greedy option, proximity cutoff 4.0, an IK-2 scalar-coupling Liouville basis at proximity level 1, and S3 symmetry groups on spins 14–16, 17–19 and 20–22. The sequence uses angle pi/2, offset 1200, sweep [2000 2000], 256 points and 512 zero-fill points on each axis. Cosine apodisation precedes the shifted 2D FFT; the plotted spectrum is the magnitude in positive mode.

## Implementation structure

The function defines the 22 proton sites and listed couplings, creates and reduces the spin system with the specified basis, and runs `liquid(...,@ct_cosy,...,'nmr')`. It then windows, Fourier-transforms and plots the result.
