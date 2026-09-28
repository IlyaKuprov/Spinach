# kernel/overloads/@ttclass/subsref.m

- Signature: `answer=subsref(ttrain,reference)`

## Purpose

Handle dot-property access, matrix-element extraction, and nested indexing for the tensor-train class.

## Parameters / inputs

- `ttrain` — tensor train object.
- `reference` — MATLAB subscript-reference structure.

## Outputs

- `answer` — the requested property, matrix element, or result of nested indexing.

## Implementation

For dot references, the supported properties are `ncores`, `ntrains`, `sizes`, `ranks`, `coeff`, `cores`, and `tolerance`; other field names raise an error. Parenthesis references require exactly two scalar indices. Logical indices are converted to numeric values, then row and column indices are validated against the tensor dimensions and converted to core indices. For each train, the function contracts the selected entries through the cores, multiplies by that train's coefficient, and sums the contributions. Advanced indexing is not implemented. Additional reference levels are applied recursively to the result.

## Source

D. Savostyanov and I. Kuprov, [`ttclass/subsref.m`](https://spindynamics.org/wiki/index.php?title=ttclass/subsref.m).
