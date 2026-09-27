# examples/kinetics/aziridine_exsy_1.m

- Signature: `aziridine_exsy_1()`

## Purpose

Simulates a NOESY/EXSY experiment on phenylaziridine with relatively slow chemical exchange. The script illustrates second-kind scalar relaxation induced by `14N`; under these conditions, first-kind scalar relaxation does not manifest. It is set to reproduce Figures 1b and 4b of the cited paper. All parameters except isotropic chemical shifts, exchange rates, and correlation times are taken from a DFT calculation; the source lists a calculation time of minutes.

## Physical / mathematical content

The spin system includes `1H` and quadrupolar `14N` nuclei. The simulation includes chemical exchange and the nitrogen quadrupolar interaction, then computes the NOESY/EXSY response.

## Reference

[doi:10.1002/ange.201410271](https://doi.org/10.1002/ange.201410271)

## Implementation structure

Defines the molecular coordinates, isotopes, magnetic field (11.75 T), and `14N` quadrupolar coupling; runs the NOESY experiment and processes the resulting signal for display.
