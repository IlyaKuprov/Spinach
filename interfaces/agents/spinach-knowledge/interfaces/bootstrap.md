# interfaces/bootstrap.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/interfaces/bootstrap.m`
- Signature: `spin_system=bootstrap(volume)`
- Total lines: 66

## Purpose

A minimal spin_system structure required to call many Spinach functions. Use this function if your system is not a spin system, and you simply want to use one of Spinach functions outside the context. Syntax: spin_system=bootstrap(volume)

## Physical / mathematical content

- This file belongs to the `interfaces` part of Spinach. Its role should be read together with nearby files in the same directory, which usually share a common physical regime or infrastructure purpose.

## Numerical / algorithmic content

- The file is built around the standard Spinach workflow: create the spin system, choose a basis or context, assemble operators/superoperators, then propagate or analyse the resulting dynamics.
- The file contains an explicit `grumble(...)` validator, which is Spinach convention for front-loading dimension, type, and regime checks before expensive linear-algebra work begins.
- The file also defines local helper function(s): `grumble()`. This usually means the public entry point is supported by tightly coupled validation or helper logic kept private to the file.

## Code-derived implementation details

### Comment-guided execution stages

- Lines 25-26: Default volume; implemented by `if ~exist('volume','var')`.
- Lines 30-31: Check consistency; implemented by `grumble(volume)`.
- Lines 33-34: Ghost spin; implemented by `sys.magnet=0; sys.isotopes={'G'}`.
- Lines 36-37: No spin interactions; implemented by `inter.zeeman.matrix=cell(1)`.
- Lines 40-41: No hygienic checks; implemented by `sys.disable={'hygiene'}`.
- Lines 43-44: Basis set; implemented by `bas.formalism='sphten-liouv'`.
- Lines 47-48: Spinach housekeeping; implemented by `if exist('volume','var')`.

### Control flow inferred from the code

- Line 26: conditional branch on `~exist('volume','var')`.
- Line 48: conditional branch on `exist('volume','var')`.

### Key state/data transformations

- Lines 27: computes `volume` using `volume='console'`.
- Lines 34: computes `sys.magnet` using `sys.magnet=0; sys.isotopes={'G'}`.
- Lines 37: computes `inter.zeeman.matrix` using `inter.zeeman.matrix=cell(1)`.
- Lines 38: computes `inter.coupling.matrix` using `inter.coupling.matrix=cell(1)`.
- Lines 41: computes `sys.disable` using `sys.disable={'hygiene'}`.
- Lines 44: computes `bas.formalism` using `bas.formalism='sphten-liouv'`.
- Lines 45: computes `bas.approximation` using `bas.approximation='none'`.
- Lines 49: computes `sys.output` using `sys.output=volume`.
- Lines 51: computes `spin_system` using `spin_system=create(sys,inter)`.

### Local helper functions

- Line 57: `grumble()` — `function grumble(volume)`. "Without love, humans would be... rare." Terry Pratchett
  - Representative operation: `if ~ischar(volume)`.
  - Representative operation: `error('volume must be a character string.')`.

## Parameters / inputs

- volume -the content of sys.output in the bootstrap
- call, set to 'hush' to make the resulting
- object suppress console output

## Outputs

- spin_system -empty but valid Spinach spin system
- description object

## Implementation structure

- A minimal spin_system structure required to call many
- Spinach functions. Use this function if your system is
- not a spin system, and you simply want to use one of
- Spinach functions outside the context. Syntax:
- spin_system=bootstrap(volume)
- volume -the content of sys.output in the bootstrap
- call, set to 'hush' to make the resulting
- object suppress console output
- spin_system -empty but valid Spinach spin system
- description object
- Default volume
- Check consistency

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `exist()`, `grumble()`, `create()`, `basis()`, `ischar()`.
