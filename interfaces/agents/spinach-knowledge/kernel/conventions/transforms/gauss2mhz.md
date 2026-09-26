# kernel/conventions/transforms/gauss2mhz.m

- Signature: `hfc_mhz=gauss2mhz(hfc_gauss,g)`

## Purpose

Converts hyperfine couplings from Gauss to MHz (linear frequency). The Gauss specification may be defined as "the magnetic field at which the electron frequency is equal to the frequency provided". Syntax: hfc_mhz=gauss2mhz(hfc_gauss,g) Arrays of any dimensions are supported. Parameters: hfc_gauss -an array of values in Gauss g -electron g-factor; if this parameter is skipped, free electron g-factor is used for conve

## Physical / mathematical content

- Convention and tensor-transform utilities. They convert among tensor parameterisations, coordinate systems, and unit systems; the underlying mathematics is linear algebra on rank-2 tensors and rotation representations.

## Numerical / algorithmic content

## Outputs

- hfc_mhz -an array of values in MHz

## Implementation structure

- Converts hyperfine couplings from Gauss to MHz (linear
- frequency). The Gauss specification may be defined as
- "the magnetic field at which the electron frequency is
- equal to the frequency provided". Syntax:
- hfc_mhz=gauss2mhz(hfc_gauss,g)
- Arrays of any dimensions are supported. Parameters:
- hfc_gauss -an array of values in Gauss
- g -electron g-factor; if this parameter
- is skipped, free electron g-factor is
- used for conversion
- hfc_mhz -an array of values in MHz
- Set the defaults
