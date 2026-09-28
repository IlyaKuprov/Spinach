# examples/nmr_liquids/inad_cyprinol.m

- Signature: `inad_cyprinol()`

## Purpose

INADEQUATE spectrum of cyprinol. The sequence selects double- quantum coherence from coupled 13C pairs and converts it back for detection. A parallel sum over isotopomers that have adjacent 13C spins is used. Calculation time: minutes

## Physical / mathematical content

The script simulates a 1D INADEQUATE spectrum of cyprinol. It generates isotopomers containing two 13C nuclei, checks the coupling between those nuclei, and only simulates pairs passing the source's `abs(J)>2*pi*1.0` threshold. The sequence selects double-quantum coherence and converts it for detection, as described in the source comment.

## Numerical / algorithmic content

The IK-1 sphten-liouv basis uses scalar-coupling connectivity, proximity level 1 and interaction level 4. The 13C acquisition uses J=50, decouples 1H, has a 10,000 Hz sweep, 5,000 Hz offset, 4,096 points and 8,192-point zero filling; an exponential apodisation parameter of 6 precedes the Fourier transform.

## Implementation structure

The field is 11.7 T. Isotopomer simulations are accumulated in a `parfor` loop; only pairs that pass the coupling threshold are built in the specified basis and passed to `liquid` with the INADEQUATE sequence. The 1D real spectrum is plotted with the configured inverted axis.
