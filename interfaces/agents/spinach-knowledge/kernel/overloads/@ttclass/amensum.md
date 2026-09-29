# kernel/overloads/@ttclass/amensum.m

Direct source: [kernel/overloads/@ttclass/amensum.m](https://github.com/IlyaKuprov/Spinach/blob/main/kernel/overloads/@ttclass/amensum.m)
Wiki: [ttclass/amensum.m](https://spindynamics.org/wiki/index.php?title=ttclass/amensum.m)

- Signature: `y=amensum(x,tol,opts)`

## Purpose and representation

Compress a sum of buffered rank-one tensor terms into one tensor train using an AMEn-style alternating projection. The input `x` and result `y` are `ttclass` objects, but they are not interchangeable representations: the input stores `N` rank-one terms in `x.cores` (`d` modes by `N` terms; each mode has its two physical dimensions and a stacked term axis), with per-mode dimensions in `x.sizes` and term coefficients in `x.coeff`; the output is one TT train whose cores have left bond rank, two physical mode axes, and right bond rank, with bond ranks that may exceed one. This is a polyadic/CP-like sum of rank-one factors being compressed to TT form, not a sum of arbitrary, already-compressed TT trains. See also [`amensolve`](amensolve.md), which requires a single TT train for its matrix and right-hand side.

## Source-defined method and checks

The method proceeds only when every entry of `x.ranks` is one. It forms each term's mode factor from `x.cores{k,i}`, contracts the buffered factors and coefficients through left/right interfaces, and updates one local TT core at a time. Optional random residual enrichment adds an auxiliary TT subspace. The largest relative local-core update controls stopping; iteration ends when it is below `tol` or `opts.max_swp` is reached. If the input has non-unit TT ranks, the current implementation raises `TT not ready yet` rather than summing general TT trains.

The source documents a target Frobenius error below `tol` times the Frobenius norm of `x`; the executable stopping test is the largest relative per-core update, so that test should not be mistaken for an independently verified global error bound.

## Options and output

`tol` is required. If `opts` is omitted, defaults are `max_swp=100`, `init_guess_rank=2`, `enrichment_rank=4`, and `verb=0`; setting `enrichment_rank=0` disables enrichment. `y` is one TT train with the same number of modes and physical mode dimensions as the input terms. Its TT bond ranks are selected during iteration; the routine records tolerance metadata on the result.
