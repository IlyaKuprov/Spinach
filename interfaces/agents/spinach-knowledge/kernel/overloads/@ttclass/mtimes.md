# kernel/overloads/@ttclass/mtimes.m

Direct mapped source: [kernel/overloads/@ttclass/mtimes.m](https://github.com/IlyaKuprov/Spinach/blob/main/kernel/overloads/@ttclass/mtimes.m) · [existing Spinach Wiki entry](https://spindynamics.org/wiki/index.php?title=ttclass/mtimes.m)

## Behaviour

The method has three operation families:

- With a `double` scalar on either side of a `ttclass`, it copies the train container, scales coefficients by the scalar, and scales tolerances by its absolute value. Cores and ranks are unchanged; this branch neither filters zero coefficients nor calls `shrink`.
- For a `ttclass` times a non-scalar `double` operand, it requires the product of the train's column-mode sizes to equal the operand's first dimension. It contracts the train cores right-to-left against that operand and accumulates each train component weighted by its coefficient. The result is a dense array with the product of the row-mode sizes as its row count and `size(b,2)` columns; it is materialised, not a tensor train. This branch has no train-rank output.
- For two `ttclass` operands, it requires equal core counts and matching contracted mode sizes (each left column-mode size equals the corresponding right row-mode size). For each pair of train components, each output core contracts the left core's right physical index with the right core's left physical index. Its dimensions are `[rA_k*rB_k, mA_k, nB_k, rA_{k+1}*rB_{k+1}]`: the uncontracted physical modes remain, and the paired TT ranks multiply. Pair coefficients multiply, producing up to `nTrainA*nTrainB` components. The method assigns each nonzero pair a tolerance based on that term's Frobenius norm and the two input relative tolerances; zero-coefficient pairs get zero tolerance. It filters zero coefficients, uses `0*unit_like(c)` if none remain, then calls `shrink(c)`, so this branch returns a recompressed tensor train whose ranks may be reduced.

The scalar branches accept only `double` scalars; a non-scalar double is supported only on the right of a tensor train, and the tensor-train pair branch is the only other supported combination. Other combinations raise an error. Core products and scalar products do not conjugate either operand; conjugate transpose is not part of this method.

No numeric example or DOI is present in the source or either existing page.
