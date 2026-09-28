# examples/kinetics/aziridine_exsy_2.m

- Signature: `aziridine_exsy_2()`

## Purpose

Simulates a NOESY/EXSY experiment on phenylaziridine in the intermediate-exchange regime, where lines broaden and first-kind scalar relaxation must be included. The script is set to reproduce Figures 1a and 4a of the cited paper. All parameters except isotropic chemical shifts, exchange rates, and correlation times are taken from a DFT calculation; the source lists a calculation time of hours, faster on a GPU.

## Physical / mathematical content

The spin system includes `1H` and quadrupolar `14N` nuclei. The model combines chemical exchange with second-kind scalar relaxation induced by nitrogen and accounts for first-kind scalar relaxation in this regime.

## Reference

[doi:10.1002/ange.201410271](https://doi.org/10.1002/ange.201410271)

## Implementation structure

Defines the molecular coordinates, isotopes, magnetic field, and `14N` quadrupolar coupling; runs the NOESY experiment and processes the resulting signal for display.
