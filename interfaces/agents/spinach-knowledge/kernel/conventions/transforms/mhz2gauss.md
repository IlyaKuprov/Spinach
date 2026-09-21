# kernel/conventions/transforms/mhz2gauss.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/kernel/conventions/transforms/mhz2gauss.m`
- Signature: `hfc_gauss=mhz2gauss(hfc_mhz,g)`
- Total lines: 58

## Purpose

Converts hyperfine couplings from MHz (linear frequency) to Gauss. The Gauss specification may be defined as "the magnetic field at which the electron frequency is equal to the frequency provided". Syntax: hfc_gauss=mhz2gauss(hfc_mhz,g) Arrays of any dimensions are supported. Parameters: hfc_mhz -an array of values in MHz g -electron g-factor; if this parameter is skipped, free electron g-factor is used for conversio

## Physical / mathematical content

- Convention and tensor-transform utilities. They convert among tensor parameterisations, coordinate systems, and unit systems; the underlying mathematics is linear algebra on rank-2 tensors and rotation representations.

## Numerical / algorithmic content

- The file contains an explicit `grumble(...)` validator, which is Spinach convention for front-loading dimension, type, and regime checks before expensive linear-algebra work begins.
- The file also defines local helper function(s): `grumble()`. This usually means the public entry point is supported by tightly coupled validation or helper logic kept private to the file.

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

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `exist()`, `grumble()`.
