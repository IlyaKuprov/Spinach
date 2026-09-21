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
