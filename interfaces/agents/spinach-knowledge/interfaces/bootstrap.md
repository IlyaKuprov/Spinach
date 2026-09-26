# interfaces/bootstrap.m

- Signature: `spin_system=bootstrap(volume)`

## Purpose

A minimal spin_system structure required to call many Spinach functions. Use this function if your system is not a spin system, and you simply want to use one of Spinach functions outside the context. Syntax: spin_system=bootstrap(volume)

## Physical / mathematical content

## Numerical / algorithmic content

- The file is built around the standard Spinach workflow: create the spin system, choose a basis or context, assemble operators/superoperators, then propagate or analyse the resulting dynamics.

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
