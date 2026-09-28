# interfaces/bootstrap.m

- Signature: `spin_system=bootstrap(volume)`

## Purpose

A minimal spin_system structure required to call many Spinach functions. Use this function if your system is not a spin system, and you simply want to use one of Spinach functions outside the context. Syntax: spin_system=bootstrap(volume)

## Physical / mathematical content

## Numerical / algorithmic content

Defaults `volume` to `'console'`, verifies that it is a character string, then constructs a ghost-spin system with no interactions and uses `create` and `basis` to return it in the requested basis.

## Parameters / inputs

- volume -the content of sys.output in the bootstrap
- call, set to 'hush' to make the resulting
- object suppress console output

## Outputs

- spin_system -empty but valid Spinach spin system
- description object

## Implementation structure

Sets `sys.magnet=0` and `sys.isotopes={'G'}`, initializes empty Zeeman and coupling matrices, disables hygiene checks, sets the output destination, then calls `create` and `basis`. A local validator rejects a non-character `volume`.
