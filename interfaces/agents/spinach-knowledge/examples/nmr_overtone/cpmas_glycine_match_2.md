# examples/nmr_overtone/cpmas_glycine_match_2.m

- Signature: `cpmas_glycine_match_2()`

## Purpose

Computes a glycine 14N-overtone/proton cross-polarisation Hartmann–Hahn profile as a function of spinning rate and 1H RF power under MAS. The source estimates a calculation time of hours and credits Ilya Kuprov, M. Carravetta, and M. Concistre.

## Physical / mathematical content

The source attributes the glycine quadrupolar tensor data to O'Dell and Ratcliffe ([DOI](http://dx.doi.org/10.1016/j.cplett.2011.08.030)). It models 14N and 1H at 14.10220742 T, with a 14N quadrupolar tensor specified by 1.18 MHz and η=0.53 and a 1H shift of 32.4. Relaxation is damped at rate 1000; the sphten-liouv basis is used without approximation. The calculation evaluates the 14N overtone signal.

## Numerical / algorithmic content

The source scans 50 proton RF powers from 10 to 200 kHz and 50 spinning rates from 20 to 90 kHz, using a 200-point octahedral powder grid and rank 5. For each pair it sets the MAS rate to the negative spinning rate, adjusts the nitrogen RF frequency to `8e3-2*rate`, and sweeps ±4 kHz around that frequency. The RF duration is 100 μs; spectra use 256 points and are reduced to summed real intensity for the 2D map.

## Implementation structure

A `parfor` loop evaluates all RF-power/spinning-rate pairs with `singlerot` and `@overtone_cp`. The resulting intensity matrix is plotted against 1H nutation frequency and sample spinning rate, both in kHz.
