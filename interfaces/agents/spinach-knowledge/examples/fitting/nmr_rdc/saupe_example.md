# examples/fitting/nmr_rdc/saupe_example.m

- Signature: `saupe_example()`

## Purpose

Extracts a Saupe order matrix from NH residual dipolar coupling (RDC) data. The experimental measurements are credited to Andras Boeszoermenyi, Thibault Viennet, and Hari Arthanari.

## Physical and mathematical content

The workflow uses NH RDC measurements and the corresponding atom coordinates from a PDB structure to fit the Saupe order matrix. It then back-calculates RDCs from the fitted matrix for comparison with the input data.

## Numerical and algorithmic content

The script reads the structure and RDC data, builds an isotope table, matches measured couplings to the relevant atom coordinates, and calls the RDC fitting routine. It plots the calculated and measured values.

## Implementation structure

The function reads the PDB and RDC inputs, maps each coupling to its pair of atoms and extracts their coordinates, fits the Saupe matrix, back-calculates the couplings, and produces the plots.
