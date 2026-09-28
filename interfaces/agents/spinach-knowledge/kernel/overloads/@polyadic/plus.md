# kernel/overloads/@polyadic/plus.m

- Signature: `c=plus(a,b)`

## Purpose

Adds polyadic or numeric operands without immediately expanding the represented Kronecker products; the result is simplified after combining terms.

## Physical / mathematical content

Addition checks matrix dimensions when both operands are non-scalar, handles zero operands as shortcuts, combines terms or wraps operands when buffers are present, then calls `simplify`.

## Numerical / algorithmic content

## Parameters / inputs

- a, b: polyadic objects or a numeric matrix

## Outputs

- c: polyadic object
- Note: use this operation sparingly -the additions are simply
- buffered, and all subsequent operations will be slower.

## Implementation structure

- Adds operands as terms in a polyadic representation without opening their Kronecker products.
- Calls `simplify` on the result.
- a,b -polyadic objects
- c -polyadic object
- Note: use this operation sparingly -the additions are simply
- buffered, and all subsequent operations will be slower.
- Check consistency
- Run shortcuts
- Possible cases
- Matrix + polyadic
