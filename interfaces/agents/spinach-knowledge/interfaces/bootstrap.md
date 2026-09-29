# interfaces/bootstrap.m

- Signature: `spin_system=bootstrap(volume)`
- Source: [interfaces/bootstrap.m](https://github.com/IlyaKuprov/Spinach/blob/main/interfaces/bootstrap.m)
- Wiki: [bootstrap.m](https://spindynamics.org/wiki/index.php?title=bootstrap.m)

## Purpose

Builds a minimal, valid Spinach system for calling Spinach functions that require a `spin_system` when no physical spin system is being modelled. It is a convenience object, not a populated molecular system.

## Input

- `volume` is optional and must be a MATLAB character array (checked with `ischar`). If omitted, it defaults to `'console'`. The value is copied to `sys.output`; the source header documents `'hush'` as the option that suppresses output.

## Construction and return value

The routine sets `sys.magnet=0` and `sys.isotopes={'G'}` (a ghost spin), and provides empty Zeeman and coupling matrices through `inter.zeeman.matrix=cell(1)` and `inter.coupling.matrix=cell(1)`. It disables the hygiene check via `sys.disable={'hygiene'}`. The basis specification is `bas.formalism='sphten-liouv'` and `bas.approximation='none'`.

It calls `create(sys,inter)`, then `basis(spin_system,bas)`, and returns that basis-initialised object as `spin_system`. The single ghost spin and absent interactions make this a structural placeholder; they do not encode a sample's spin Hamiltonian.

## Dependencies and guardrails

The implementation depends on Spinach's `create` and `basis` functions. Its local validator rejects non-character `volume` values; it does not accept a MATLAB string scalar under that check.
