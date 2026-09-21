# kernel/utilities/chemshifts.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/kernel/utilities/chemshifts.m`
- Signature: `[cs_ppm,cs_hz]=chemshifts(spin_system)`
- Total lines: 60

## Purpose

Returns the chemical shifts of every spin in the system relative to the carrier frequency in the current magnet. Syntax: [cs_ppm,cs_hz]=chemshifts(spin_system)

## Physical / mathematical content

- General mathematical and infrastructure utilities. This area contains finite differences, perturbation theory, graph algorithms, spectral densities, tensor algebra, hash/report helpers, and other reusable numerical components.

## Numerical / algorithmic content

- The file contains an explicit `grumble(...)` validator, which is Spinach convention for front-loading dimension, type, and regime checks before expensive linear-algebra work begins.
- The file also defines local helper function(s): `grumble()`. This usually means the public entry point is supported by tightly coupled validation or helper logic kept private to the file.

## Parameters / inputs

- spin_system -spin system descriptor object

## Outputs

- cs_ppm -chemical shifts in ppm
- cs_hz -chemical shifts in Hz

## Implementation structure

- Returns the chemical shifts of every spin in the system relative
- to the carrier frequency in the current magnet. Syntax:
- [cs_ppm,cs_hz]=chemshifts(spin_system)
- spin_system -spin system descriptor object
- cs_ppm -chemical shifts in ppm
- cs_hz -chemical shifts in Hz
- Check consistency
- Preallocate outputs
- Fill the outputs
- Get isotropic Zeeman frequencies
- Subtract carrier
- Convert into ppm

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `grumble()`, `cs_ppm()`, `cs_hz()`, `isfield()`.
