# kernel/overloads/@polyadic/inflate.m

MATLAB source: <https://github.com/IlyaKuprov/Spinach/blob/main/kernel/overloads/@polyadic/inflate.m>
Spinach Wiki: <https://spindynamics.org/wiki/index.php?title=polyadic/inflate.m>

## Meaning

A polyadic value stores a sum of terms in `p.cores`: each `p.cores{n}` is one term, and its cells are the ordered matrix factors in that term's Kronecker product. The `prefix` and `suffix` cell arrays hold matrix factors acting on the left and right of the core sum. This routine converts that representation to an explicit matrix; it is the materialisation counterpart to building or composing the factorised representation.

## Construction and dimensions

Nested polyadics in core factors, prefixes, and suffixes are first recursively passed to `inflate`. The output row and column counts for the core sum are products of the row and column sizes of the factors in the first term, respectively. For each term, the function forms a left-to-right Kronecker product of its factors, extracts its nonzero row indices, column indices, and values with `find`, combines those triplets across terms, and constructs a sparse matrix of the computed dimensions. If there are no triplets, it constructs an empty sparse matrix with those dimensions.

The dimensions are taken from `p.cores{1}`; this file does not add an explicit consistency check for the dimensions of later terms. The source does not implement broadcasting or dimension repair: products and affix actions use MATLAB's matrix operations and therefore require compatible dimensions.

## Operator composition and storage

After assembling the core sum, prefixes are left-multiplied in reverse cell order, so the resulting product is `prefix{1} * prefix{2} * core_sum` when there are two prefixes. Suffixes are right-multiplied in forward cell order, giving `core_sum * suffix{1} * suffix{2}`. The core accumulator is explicitly made sparse from triplets. The source header notes that full prefix or suffix factors may produce a full result; the function does not convert the final product back to sparse storage.

No orientation variable or orientation-specific action is present in this function: it composes the matrix factors stored in the polyadic value.

Function handle cores are rejected before expansion, including inside nested polyadics; use operator actions rather than materialising their matrices. A singleton affix can produce an `opium` identity during multiplication; the routine converts that result to an explicit sparse matrix before returning.
