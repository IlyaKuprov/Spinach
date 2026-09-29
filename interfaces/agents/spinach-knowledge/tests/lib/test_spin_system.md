# tests/lib/test_spin_system.m

## Purpose

Builds a small quiet Spinach spin system for tests.

## Behaviour

The function applies quiet settings used by regression tests before constructing the spin system:

- Sets `sys.output` to `'hush'`.
- Appends `'hygiene'` to `sys.disable` if that field already exists (using `unique` on the concatenated cell array), otherwise initialises `sys.disable` to `{'hygiene'}`.
- Sets `sys.parallel` to `{'local',1}`.
- Sets `sys.parprops` to `{}`.

It then builds the Spinach object with `create(sys,inter)` and the basis with `basis(spin_system,bas)`, returning the resulting spin system object.

## Inputs and outputs

Syntax:

```
spin_system=test_spin_system(sys,inter,bas)
```

Inputs:

- `sys` — Spinach system specification.
- `inter` — Spinach interaction specification.
- `bas` — Spinach basis specification.

Output:

- `spin_system` — Spinach spin system object.

## References

- [Source file on GitHub](https://github.com/IlyaKuprov/Spinach/blob/main/tests/lib/test_spin_system.m)
