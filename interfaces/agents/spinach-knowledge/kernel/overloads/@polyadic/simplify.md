# kernel/overloads/@polyadic/simplify.m

- Signature: `p=simplify(p)`

## Purpose

Recursively simplifies prefix and suffix buffers and core terms: it removes zero and identity factors, absorbs scalar or `opium` factors, flattens eligible single-core nested polyadics, and combines adjacent `opium` factors. A zero representation becomes a sparse zero matrix; a bare single-core representation is returned directly.

## Physical / mathematical content

Simplification preserves the matrix represented by the polyadic while reducing redundant or zero structure.

## Numerical / algorithmic content

## Parameters / inputs

- `p`: a polyadic object

## Outputs

- `p`: a polyadic or numeric object

## Implementation structure

- Validates the input and obtains its represented dimensions.
- Repeats simplification until no changes remain, processing prefixes, suffixes, and core terms recursively.
- Returns sparse zero output for an empty/zero representation, and unwraps a bare single-core polyadic.
