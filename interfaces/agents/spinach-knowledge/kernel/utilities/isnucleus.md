# kernel/utilities/isnucleus.m

## Purpose

Returns `true` if a given spin specification string is a nucleus, and `false` otherwise. The function is a lightweight classifier used to distinguish nuclear spin specifications from electron-related and other non-nuclear entries in Spinach spin system specifications.

Source: <https://github.com/IlyaKuprov/Spinach/blob/main/kernel/utilities/isnucleus.m>

## Behaviour

1. The input is first passed to an internal consistency-checking subfunction `grumble`, which:
   - Errors with `'spin_spec must be a character string.'` if the input is not a character string (`~ischar(spin_spec)`).
   - Calls `spin(spin_spec)` to verify that the specification is a valid spin specification; an invalid specification causes `spin` to throw its own error.
2. After validation, classification is a simple name-matching check:
   - If the first character of `spin_spec` is one of `'E'`, `'C'`, `'V'`, or `'T'`, the verdict is `false`.
   - If the whole string equals `'G'`, `'E'`, `'N'`, or `'M'`, the verdict is `false`.
   - Otherwise, the verdict is `true`.
3. The returned values are produced with `false()` and `true()`, so the output is a logical scalar.

The single-letter prefixes handled as non-nuclear correspond to Spinach's electron and pseudo-spin specification families (including the `E` electron prefix); the check is purely string-based and involves no database lookup beyond the validity check performed by `spin`.

## Inputs and outputs

**Syntax**

```matlab
verdict = isnucleus(spin_spec)
```

**Parameters**

- `spin_spec` — a character string containing a spin specification. Must be a valid specification accepted by `spin`.

**Outputs**

- `verdict` — logical scalar; `true` for a nucleus, `false` otherwise.

**Errors**

- `'spin_spec must be a character string.'` — raised when the input is not a character string.
- Errors from `spin(spin_spec)` — raised when the specification string is not a valid spin specification.

## References

- Spinach Wiki page for this function: <https://spindynamics.org/wiki/index.php?title=isnucleus.m>
- Source file on GitHub: <https://github.com/IlyaKuprov/Spinach/blob/main/kernel/utilities/isnucleus.m>
