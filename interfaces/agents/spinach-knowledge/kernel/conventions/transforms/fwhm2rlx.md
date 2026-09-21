# kernel/conventions/transforms/fwhm2rlx.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/kernel/conventions/transforms/fwhm2rlx.m`
- Signature: `r2rate=fwhm2rlx(fwhm)`
- Total lines: 44

## Purpose

Converts full width at half-maximum (FWHM) of an NMR signal into an approximation of the R2 rate. Syntax: r2rate=fwhm2rlx(fwhm)

## Physical / mathematical content

- Convention and tensor-transform utilities. They convert among tensor parameterisations, coordinate systems, and unit systems; the underlying mathematics is linear algebra on rank-2 tensors and rotation representations.

## Numerical / algorithmic content

- The file contains an explicit `grumble(...)` validator, which is Spinach convention for front-loading dimension, type, and regime checks before expensive linear-algebra work begins.
- The file also defines local helper function(s): `grumble()`. This usually means the public entry point is supported by tightly coupled validation or helper logic kept private to the file.

## Parameters / inputs

- fwhm -full width at half-maximum, Hz

## Outputs

- r2rate -approximate R2 relaxation rate, Hz
- Note: FWHM is not a reliable measure of the transverse
- relaxation rate. The value obtained from this
- function should be treated as an upper bound.
- Note: Lorentzian line shape is assumed.

## Implementation structure

- Converts full width at half-maximum (FWHM) of an NMR
- signal into an approximation of the R2 rate. Syntax:
- r2rate=fwhm2rlx(fwhm)
- fwhm -full width at half-maximum, Hz
- r2rate -approximate R2 relaxation rate, Hz
- Note: FWHM is not a reliable measure of the transverse
- relaxation rate. The value obtained from this
- function should be treated as an upper bound.
- Note: Lorentzian line shape is assumed.
- Check consistency
- Run the conversion
- Consistency enforcement

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `grumble()`, `any()`, `fwhm()`.
