# kernel/conventions/transforms/mt2hz.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/kernel/conventions/transforms/mt2hz.m`
- Signature: `hfc_hz=mt2hz(hfc_mt,g)`
- Total lines: 59

## Purpose

Converts hyperfine couplings from milliTesla to Hz (linear frequency). The milliTesla specification may be defined as "the magnetic field at which the electron frequency is equ- al to the frequency provided". Syntax: hfc_hz=mt2hz(hfc_mt,g) Arrays of any dimension are supported. Parameters: hfc_mt -an array of values in mT g -electron g-factor; if this parameter is skipped, free electron g-factor is used for conversio

## Physical / mathematical content

- Convention and tensor-transform utilities. They convert among tensor parameterisations, coordinate systems, and unit systems; the underlying mathematics is linear algebra on rank-2 tensors and rotation representations.

## Numerical / algorithmic content

- The file contains an explicit `grumble(...)` validator, which is Spinach convention for front-loading dimension, type, and regime checks before expensive linear-algebra work begins.
- The file also defines local helper function(s): `grumble()`. This usually means the public entry point is supported by tightly coupled validation or helper logic kept private to the file.

## Outputs

- hfc_hz -an array of values in Hz

## Implementation structure

- Converts hyperfine couplings from milliTesla to Hz (linear
- frequency). The milliTesla specification may be defined as
- "the magnetic field at which the electron frequency is equ-
- al to the frequency provided". Syntax:
- hfc_hz=mt2hz(hfc_mt,g)
- Arrays of any dimension are supported. Parameters:
- hfc_mt -an array of values in mT
- g -electron g-factor; if this parameter
- is skipped, free electron g-factor is
- used for conversion
- hfc_hz -an array of values in Hz
- Set the defaults

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `exist()`, `grumble()`.
