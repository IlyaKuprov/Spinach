# kernel/overloads/@ttclass/full.m

- Signature: `answer=full(ttrain)`

## Purpose

Materialises a tensor-train matrix as a dense numeric matrix. A TT stores each sum term as a column of core cells, with a coefficient for that term; each core carries left/right bond ranks and one row-mode and one column-mode dimension. This is a tensor-train representation, not a polyadic-object storage interface.

## Contraction and result shape

The result is initialised with `zeros(size(ttrain))`, so its matrix shape is the product of the row-mode sizes by the product of the column-mode sizes. For each train term, the function starts from the last core, contracts toward the first using the adjacent bond ranks, and reshapes/permutes the accumulated row and column mode axes into matrix order. It multiplies that dense term by its stored coefficient and adds it to the result.

The allocation is dense and can be very large. The function contains no explicit class, rank-consistency, or mode-size validation; it relies on the object's core, rank, size, and `size` methods to describe a valid TT.

## Sources

- [MATLAB source](https://github.com/IlyaKuprov/Spinach/blob/main/kernel/overloads/@ttclass/full.m)
- [Spinach Wiki](https://spindynamics.org/wiki/index.php?title=ttclass/full.m)
