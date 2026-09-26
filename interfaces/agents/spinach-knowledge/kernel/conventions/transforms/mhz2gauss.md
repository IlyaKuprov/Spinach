# kernel/conventions/transforms/mhz2gauss.m

- Signature: `hfc_gauss=mhz2gauss(hfc_mhz,g)`

## Purpose

Converts hyperfine couplings from MHz (linear frequency) to Gauss. The Gauss specification may be defined as "the magnetic field at which the electron frequency is equal to the frequency provided". Syntax: hfc_gauss=mhz2gauss(hfc_mhz,g) Arrays of any dimensions are supported. Parameters: hfc_mhz -an array of values in MHz g -electron g-factor; if this parameter is skipped, free electron g-factor is used for conversio

## Physical / mathematical content

- Convention and tensor-transform utilities. They convert among tensor parameterisations, coordinate systems, and unit systems; the underlying mathematics is linear algebra on rank-2 tensors and rotation representations.

## Numerical / algorithmic content

## Outputs

- hfc_gauss -an array of values in Gauss

## Implementation structure

- Converts hyperfine couplings from MHz (linear frequency)
- to Gauss. The Gauss specification may be defined as
- "the magnetic field at which the electron frequency is
- equal to the frequency provided". Syntax:
- hfc_gauss=mhz2gauss(hfc_mhz,g)
- Arrays of any dimensions are supported. Parameters:
- hfc_mhz -an array of values in MHz
- g -electron g-factor; if this parameter
- is skipped, free electron g-factor is
- used for conversion
- hfc_gauss -an array of values in Gauss
- Set the defaults
