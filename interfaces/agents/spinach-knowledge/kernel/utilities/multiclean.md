# kernel/utilities/multiclean.m

## Purpose

`C=multiclean(spin_system,A,B,tol)` multiplies matrices in output-column panels and cleans each panel before assembly, reducing transient fill-in. Complete inner products are retained. CPU and GPU, full and sparse inputs are supported; eligible sparse GPU products use the CUDA backend internally.

## Numerical contract

`tol` must be a finite, real, non-negative floating-point scalar. NaN, infinity and integer tolerances are errors. Zero tolerance and the cleanup-disable switch are handled by `clean_up`, without separate bypass logic in the helper. The supplied spin system is neither copied nor modified.

Panel nonzero entries are assembled into a sparse matrix and passed through ordinary `clean_up` for the final storage choice. Repeated rounding and blocked multiplication can differ from one native product by numerical roundoff; bitwise preservation is not promised. Scalars, empty products and symbolic operator objects use native multiplication followed by cleanup.

Panel width is based on input nonzero counts and the storage for three dense output panels. This favours sparse-workload memory over speed; it is not a bound on backend workspace, final output or concatenation storage.
