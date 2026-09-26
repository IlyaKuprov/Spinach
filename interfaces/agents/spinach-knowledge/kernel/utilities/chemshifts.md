# kernel/utilities/chemshifts.m

- Signature: `[cs_ppm,cs_hz]=chemshifts(spin_system)`

## Purpose

Returns the chemical shifts of every spin in the system relative to the carrier frequency in the current magnet. Syntax: [cs_ppm,cs_hz]=chemshifts(spin_system)

## Physical / mathematical content

- General mathematical and infrastructure utilities. This area contains finite differences, perturbation theory, graph algorithms, spectral densities, tensor algebra, hash/report helpers, and other reusable numerical components.

## Numerical / algorithmic content

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
