# examples/nmr_overtone/cpmas_glycine_match_1.m

- Signature: `cpmas_glycine_match_1()`

## Purpose

Maps the glycine 14N-overtone/proton cross-polarisation Hartmann–Hahn profile against 1H RF power under MAS, using a rough powder grid. The source estimates minutes of computation and credits Ilya Kuprov, M. Carravetta, and M. Concistre.

## Physical / mathematical content

The source attributes the glycine quadrupolar tensor data to O'Dell and Ratcliffe ([DOI](http://dx.doi.org/10.1016/j.cplett.2011.08.030)). The model uses 14N and 1H at 14.10220742 T, quadrupolar parameters 1.18 MHz and η=0.53 for 14N, a 1H shift of 32.4, damping rate 300, and the sphten-liouv basis without approximation. The MAS rate is −19.840 kHz and the nitrogen overtone RF frequency is 48 kHz.

## Numerical / algorithmic content

The rough powder grid is `rep_2ang_200pts_oct`; spectra use rank 7, 256 points over [44, 52] kHz. The 1H RF power is sampled at 15 values from 25 to 39 kHz, while the 14N RF power is 55 kHz; the RF duration is 100 μs. Each setting is simulated with `singlerot` and `@overtone_cp`.

## Implementation structure

The function prepares the system, basis and CP operators, loops over the 15 proton-power settings, computes a spectrum for each, and displays the spectra in a row of panels.
