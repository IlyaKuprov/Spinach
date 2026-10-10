# kernel/overloads/@ttclass/amensolve.m

Direct source: [kernel/overloads/@ttclass/amensolve.m](https://github.com/IlyaKuprov/Spinach/blob/main/kernel/overloads/@ttclass/amensolve.m)
Wiki: [ttclass/amensolve.m](https://spindynamics.org/wiki/index.php?title=ttclass/amensolve.m)

- Signature: `x=amensolve(A,y,tol,opts,x0)`

## Purpose

Approximately solve `A*x=y` by alternating minimal energy (AMEn) sweeps on tensor-train (TT) cores.

## TT representation, dimensions, and checks

`A`, `y`, and `x0` are `ttclass` objects. `A` stores matrix modes: each row of `A.sizes` holds that core's row-mode and column-mode dimensions. The source requires equal row- and column-mode sizes at every core. It requires `y` to have the same number of cores, row-mode sizes matching the matrix column-mode sizes, and singleton second mode at every core, so `y` is a vector TT. If supplied, `x0` must also be a `ttclass` vector TT with matching core count and mode sizes. The solver requires `A` and `y` to be shrunk to one train (`ntrains <= 1`); it is not a solver for a buffered/polyadic sum of trains. It also checks that `tol` is the relative approximation/stopping tolerance (the source suggests `1e-6` as a starting value); `tol` and `opts.sv_floor` are checked as finite, non-negative real scalars.

Each matrix core is stored with left bond rank, row-mode, column-mode, and right bond rank axes; a vector core has left bond rank, physical mode, and right bond rank axes. The implementation contracts left and right environment reductions across the TT bonds with a matrix core to form the local operator, solves for the local vector core, then moves the active site through the train. SVD truncation controls the solution ranks; residual enrichment can augment the train. If `x0` is omitted, a random TT initial guess is made using `opts.init_guess_rank` and orthogonalised. A zero right-hand side returns the zero TT immediately.

## Options

If `opts` is omitted or empty, fields default to: `nswp=20`, `init_guess_rank=2`, `enrichment_rank=4`, `resid_damp=2`, `rmax=Inf`, `max_full_size=500`, `local_iters=100`, `sv_floor=1e-10`, and `verb=1`. `enrichment_rank=0` disables enrichment. The direct solve is selected when the local vector block size `rx(k)*sz(k)*rx(k+1)` is below `max_full_size`; otherwise the local system uses BiCGSTAB, capped by `local_iters`. The singular-value floor and `rmax` bound SVD truncation. Sweeps continue until the measured maximum local error is below `tol` or `nswp` is reached (with the implementation completing its sweep-direction pass).

## Output

`x` is a single vector TT with the same per-mode vector dimensions as `y`; its ranks are determined by the solve and truncation. The source describes `tol` as a relative approximation and stopping tolerance; the in-code sweep stopping test tracks the maximum local error and reports a residual separately.
