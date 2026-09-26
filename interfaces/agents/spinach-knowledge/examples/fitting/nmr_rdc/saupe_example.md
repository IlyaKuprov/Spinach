# examples/fitting/nmr_rdc/saupe_example.m

- Signature: `saupe_example()`

## Purpose

Extracting Saupe order matrix from NH RDC data. Experimental measurements kindly provided by Andras Boeszoermenyi, Thibault Viennet, and Hari Arthanari.

## Physical / mathematical content

- Fitting examples. These files formulate parameter-estimation workflows in which simulated spectra or observables are matched to data, usually through nonlinear optimisation, residual construction, and physically constrained parameterisations.

## Numerical / algorithmic content

## Implementation structure

- Extracting Saupe order matrix from NH RDC data. Experimental
- measurements kindly provided by Andras Boeszoermenyi, Thibault
- Viennet, and Hari Arthanari.
- Read the PDB file
- Read RDC data
- Make isotope table
- Match up RDCs with coordinates
- Locate both atoms
- Extract coordinates
- Call RDC fitter
- Back-calculate RDCs
- Do the plotting
