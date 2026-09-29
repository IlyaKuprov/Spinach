# kernel/overloads/@ttclass/truncate.m

[MATLAB source](https://github.com/IlyaKuprov/Spinach/blob/main/kernel/overloads/@ttclass/truncate.m) · [Spin Dynamics Wiki](https://spindynamics.org/wiki/index.php?title=ttclass/truncate.m)

## Signature

`ttout=truncate(tt)`

## Purpose and input

This internal stage recompresses one tensor train by a right-to-left SVD sweep; use [shrink](shrink.md) rather than calling it directly. The input must contain one train (`tt.ntrains=1`) and be orthogonalised left-to-right. The output has that train orthogonalised right-to-left.

## Core indexing and truncation

For `d=tt.ncores`, core `k` has shape `[r(k),m(k,1),m(k,2),r(k+1)]`. The sweep visits `k=d:-1:2`. It reshapes core `k` to `[r(k),m(k,1)*m(k,2)*r(k+1)]` and core `k-1` to `[r(k-1)*m(k-1,1)*m(k-1,2),r(k)]`, takes an economy SVD of the former, and contracts the retained left singular vectors and singular values into core `k-1`. The updated bond rank is `rnew`; core `k` is rebuilt as `[rnew,m(k,1),m(k,2),r(k+1)]` and core `k-1` as `[r(k-1),m(k-1,1),m(k-1,2),rnew]`.

The requested tolerance is absolute in Frobenius norm. At each of the `d-1` cuts the routine sets `eta = tt.tolerance(1,1) / abs(tt.coeff(1,1)) / sqrt(d)` and passes `eta * norm(S,'fro')` to [frob_chop](../../utilities/frob_chop.md). That helper first zeros singular values strictly below `numel(s)*eps*max(abs(s))`; it then accumulates squared values from the smallest upward and chooses the smallest leading prefix whose discarded tail is below the requested cutoff. Equivalently, if the first tail sum reaching the cutoff is at index `k`, the retained rank is `numel(s)-k+1`. This gives the per-cut relative threshold after coefficient scaling and equal quadrature allocation across the `d` cores.

After the sweep, the first core is divided by its 2-norm and that norm is multiplied into `ttout.coeff(1,1)`. The result remains a single buffered train, with updated bond ranks and right-to-left orthogonalised cores.

## Output

- `ttout` — the recompressed tensor train.
