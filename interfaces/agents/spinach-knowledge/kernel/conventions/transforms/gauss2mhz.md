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

## Code-derived implementation details

### Comment-guided execution stages

- Lines 26-27: Set the defaults; implemented by `if ~exist('g','var')`.
- Lines 32-33: Check consistency; implemented by `grumble(hfc_gauss,g)`.
- Lines 35-36: Do the conversion; implemented by `muB=9.274009994*10^-24`.

### Control flow inferred from the code

- Line 27: conditional branch on `~exist('g','var')`.

### Key state/data transformations

- Lines 29: computes `g` using `g=2.0023193043622`.
- Lines 36: computes `muB` using `muB=9.274009994*10^-24`.
- Lines 37: computes `hbar` using `hbar=1.054571628e-34`.
- Lines 38: computes `C` using `C=1e-10*g*muB/(hbar*2*pi)`.
- Lines 39: computes `hfc_mhz` using `hfc_mhz=C*hfc_gauss`.

### Local helper functions

- Line 44: `grumble()` — `function grumble(hfc_gauss,g)`. A celibate clergy is an especially good idea, because it tends to suppress any hereditary propensity toward
  - Representative operation: `if (~isnumeric(hfc_gauss))||(~isreal(hfc_gauss))`.
  - Representative operation: `error('hfc_gauss must be an array of real numbers.')`.

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
