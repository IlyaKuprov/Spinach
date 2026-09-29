# kernel/utilities/polinfo.m

## Purpose

Prints an ASCII diagram of a `polyadic` object to the console, showing its size, prefix, Kronecker-core terms, and suffix, including nested polyadic elements.

Source: <https://github.com/IlyaKuprov/Spinach/blob/main/kernel/utilities/polinfo.m>

## Behaviour

- Signature: `polinfo(p,level,label)`, with `level` defaulting to `0` and `label` defaulting to `'polyadic'` when not supplied.
- Validates inputs via the internal `grumble` function: `p` must be a `polyadic` object, `level` a non-negative real integer scalar, and `label` a character string; violations raise errors (`'p must be a polyadic object.'`, `'level must be a non-negative integer.'`, `'label must be a character string.'`).
- Indentation is `4*level` spaces, so nested calls print deeper levels.
- Prints the label and the object's dimensions as `label [nrowsxncols]`.
- For `polyadic` inputs, prints a `prefix:` header when `p.prefix` is non-empty, then lists each prefix element: nested `polyadic` elements recurse via `polinfo(p.prefix{n},level+2,sprintf('polyad %d',n))`; `opium` objects print as `opium  n [nrowsxncols]`; anything else prints as `matrix  n [nrowsxncols]`.
- Iterates `p.cores`, printing `kron n` for each core cell and then each factor within it, with the same polyadic/opium/matrix classification and recursion at `level+2`.
- For `polyadic` inputs, prints a `suffix:` header when `p.suffix` is non-empty, then lists suffix elements with the same classification and recursion.
- The header comment notes that polyadic objects can be huge and that the code avoids making memory copies.

## Inputs and outputs

Inputs:

- `p` — polyadic object to diagram.
- `level` — optional non-negative integer indentation level (default `0`).
- `label` — optional character string label for the printed header (default `'polyadic'`).

Outputs:

- An ASCII diagram printed to the console; no return values.

## References

- Spinach Wiki: <https://spindynamics.org/wiki/index.php?title=polyadic/polinfo.m>
- Source file: <https://github.com/IlyaKuprov/Spinach/blob/main/kernel/utilities/polinfo.m>
