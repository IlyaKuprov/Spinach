# kernel/utilities/multiclean.m

## Purpose

`C=multiclean(spin_system,A,B,tol)` computes a cleaned-up matrix product. When the full product fits into the memory share of the current worker, it is ordinary `clean_up(A*B,tol)`. Otherwise the product is formed in column panels that are cleaned as they are made, which limits transient fill-in of sparse products at the cost of runtime.

## Numerical contract

`tol` must be a finite, real, non-negative floating-point scalar; NaN, infinity and integer tolerances are errors. Zero tolerance and the cleanup-disable switch are handled by `clean_up`. The supplied spin system is neither copied nor modified.

The memory share is host free memory divided by `spin_system.sys.parallel{2}`, on the assumption that all workers perform the same operation. The full product is deemed to fit when `rows*cols` elements, with two `clean_up` temporaries, are within that share. Panel width follows the same bound. GPU arrays and symbolic operator objects always use ordinary multiplication; sparse double GPU products use `cuda_sparse_by_sparse`. Panel results may differ from one native product by numerical roundoff.
