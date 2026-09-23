# kernel/utilities/min_int_type.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/kernel/utilities/min_int_type.m`
- Signature: `type=min_int_type(max_val,issigned)`
- Total lines: 106

## Purpose

Minimum integer data type sufficient to store the specified value. Useful in many indexing operations in the Spinach kernel where double precision would be a massive overkill.

## Physical / mathematical content

- General mathematical and infrastructure utilities. This area contains finite differences, perturbation theory, graph algorithms, spectral densities, tensor algebra, hash/report helpers, and other reusable numerical components.

## Numerical / algorithmic content

- The signed ladder `int8`, `int16`, `int32`, `int64` or the unsigned ladder `uint8`, `uint16`, `uint32`, `uint64` is walked upwards and the first class whose `intmax` is at least `max_val` is returned as a class name string; a value beyond the 64-bit class is an error.
- The bound is inclusive and the smallest accepted maximum is zero, which still needs one byte: `basis.m` asks for the class of the largest single-spin state index, `max(mults)^2-1`, and that is zero for a system of multiplicity-one ghost spins.
- The grumbler requires a real non-negative integer scalar and `issigned` equal to `'signed'` or `'unsigned'`.

## Syntax

```matlab
type=min_int_type(max_val,issigned)
```

## Parameters / inputs

- max_val - maximum value that the integer must cover
- issigned - whether the integer needs to cover the negative values: 'signed' or 'unsigned'

## Outputs

- type - Matlab data type to use
