# kernel/utilities/iselectron.m

## Purpose

Returns `true` if the particle specified is an electron, and `false` otherwise.

## Behaviour

- The function first validates the input via an internal consistency check (`grumble`):
  - Errors with `'spin_spec must be a character string.'` if the input is not a character string.
  - Calls `spin(spin_spec)` to verify that the specification is a valid Spinach particle specification.
- After validation, the function performs a simple matching check: if the first element of `spin_spec` is `'E'`, the verdict is `true`; otherwise it is `false`.

## Inputs and outputs

**Syntax**

```matlab
verdict = iselectron(spin_spec)
```

**Inputs**

- `spin_spec` — a Spinach particle specification (character string).

**Outputs**

- `verdict` — `true` for an electron, `false` otherwise.

## References

- Source: [kernel/utilities/iselectron.m](https://github.com/IlyaKuprov/Spinach/blob/main/kernel/utilities/iselectron.m)
- Wiki: <https://spindynamics.org/wiki/index.php?title=iselectron.m>
