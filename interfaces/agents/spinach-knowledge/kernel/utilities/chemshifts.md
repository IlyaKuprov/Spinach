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

## Code-derived implementation details

### Comment-guided execution stages

- Lines 22-23: Check consistency; implemented by `grumble(spin_system)`.
- Lines 25-26: Preallocate outputs; implemented by `cs_ppm=zeros(spin_system.comp.nspins,1)`.
- Lines 29-30: Fill the outputs; implemented by `for n=1:spin_system.comp.nspins`.
- Lines 32-33: Get isotropic Zeeman frequencies; implemented by `iso=trace(spin_system.inter.zeeman.matrix{n})/3`.
- Lines 35-36: Subtract carrier; implemented by `iso=iso-spin_system.inter.basefrqs(n)`.
- Lines 38-39: Convert into ppm; implemented by `cs_ppm(n)=1e6*iso/spin_system.inter.basefrqs(n)`.
- Lines 41-42: Convert into Hz; implemented by `cs_hz(n)=-iso/(2*pi)`.

### Control flow inferred from the code

- Line 30: `for` loop over `n=1:spin_system.comp.nspins`.

### Key state/data transformations

- Lines 26: computes `cs_ppm` using `cs_ppm=zeros(spin_system.comp.nspins,1)`.
- Lines 27: computes `cs_hz` using `cs_hz=zeros(spin_system.comp.nspins,1)`.
- Lines 33: computes `iso` using `iso=trace(spin_system.inter.zeeman.matrix{n})/3`.
- Lines 39: computes `cs_ppm(n)` using `cs_ppm(n)=1e6*iso/spin_system.inter.basefrqs(n)`.
- Lines 42: computes `cs_hz(n)` using `cs_hz(n)=-iso/(2*pi)`.

### Local helper functions

- Line 49: `grumble()` — `function grumble(spin_system)`. I have never understood why it is 'greed' to want to keep the money you have earned, but not greed to want
  - Representative operation: `if (~isfield(spin_system,'comp'))||(~isfield(spin_system,'inter'))`.
  - Representative operation: `error('the spin system object does not contain the required information.')`.

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
