# kernel/conventions/transforms/gauss2mhz.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/kernel/conventions/transforms/gauss2mhz.m`
- Signature: `hfc_mhz=gauss2mhz(hfc_gauss,g)`
- Total lines: 58

## Purpose

Converts hyperfine couplings from Gauss to MHz (linear frequency). The Gauss specification may be defined as "the magnetic field at which the electron frequency is equal to the frequency provided". Syntax: hfc_mhz=gauss2mhz(hfc_gauss,g) Arrays of any dimensions are supported. Parameters: hfc_gauss -an array of values in Gauss g -electron g-factor; if this parameter is skipped, free electron g-factor is used for conve

## Physical / mathematical content

- Convention and tensor-transform utilities. They convert among tensor parameterisations, coordinate systems, and unit systems; the underlying mathematics is linear algebra on rank-2 tensors and rotation representations.

## Numerical / algorithmic content

- The file contains an explicit `grumble(...)` validator, which is Spinach convention for front-loading dimension, type, and regime checks before expensive linear-algebra work begins.
- The file also defines local helper function(s): `grumble()`. This usually means the public entry point is supported by tightly coupled validation or helper logic kept private to the file.

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

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `exist()`, `grumble()`.
