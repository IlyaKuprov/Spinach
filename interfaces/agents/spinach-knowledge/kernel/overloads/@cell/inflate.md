# kernel/overloads/@cell/inflate.m

- Signature: `A=inflate(A)`

The cell overload calls `inflate` on each value in linear-index order and stores each returned value in the same cell. It therefore preserves the cell-array shape while allowing dispatch to materialise individual entries; it does not broadcast across cells or impose a cell-content type check.

For a polyadic entry, `@polyadic/inflate` recursively expands nested polyadics in its cores, prefixes, and suffixes. It forms each core term from the Kronecker product of that term's factors, collects the terms into a sparse matrix, then left-multiplies prefixes and right-multiplies suffixes in their stored order. This is materialisation plus ordered operator composition, not a tensor contraction or elementwise broadcast. The core matrix is assembled sparse; full prefixes or suffixes can yield a full final result. If a prefix or suffix is an `opium`, its product dispatches to `@opium/mtimes`: compatible matrix products are dimension-checked and yield a scaled numeric matrix, while an opium-opium product combines coefficients; a numeric scalar on the left scales the opium coefficient. These operations do not add a broadcasting rule.

The cell wrapper itself validates no input domain; its documented use is a cell array of polyadics. Polyadic materialisation derives the core matrix dimensions from the first core term, so the source does not provide a general shape-validation guarantee for inconsistent terms.

## References

- MATLAB source: [`kernel/overloads/@cell/inflate.m`](https://github.com/IlyaKuprov/Spinach/blob/main/kernel/overloads/@cell/inflate.m)
- Polyadic implementation: [`kernel/overloads/@polyadic/inflate.m`](https://github.com/IlyaKuprov/Spinach/blob/main/kernel/overloads/@polyadic/inflate.m)
- Operator multiplication: [`kernel/overloads/@opium/mtimes.m`](https://github.com/IlyaKuprov/Spinach/blob/main/kernel/overloads/@opium/mtimes.m)
- Spinach Wiki: [`cell/inflate.m`](https://spindynamics.org/wiki/index.php?title=cell/inflate.m); [`polyadic/inflate.m`](https://spindynamics.org/wiki/index.php?title=polyadic/inflate.m)
