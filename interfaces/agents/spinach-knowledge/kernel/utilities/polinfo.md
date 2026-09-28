# kernel/utilities/polinfo.m

- Signature: `polinfo(p,level,label)`

## Purpose

Prints an ASCII description of a polyadic object to the console.

## Parameters / inputs

- `p` — polyadic object to describe.
- `level` — optional non-negative integer scalar indentation level; defaults to `0`.
- `label` — optional character-string label; defaults to `'polyadic'`.

## Output

There is no return value; the function prints the description. It indents by four spaces per level and shows the label and object dimensions.

## Implementation structure

The description traverses the object's prefix, Kronecker terms and core entries, then its suffix. Nested polyadic entries are described recursively at `level+2`; other entries are identified by type and dimensions. The implementation avoids copying polyadic objects, which can be very large.

## Reference

- <https://spindynamics.org/wiki/index.php?title=polyadic/polinfo.m>
