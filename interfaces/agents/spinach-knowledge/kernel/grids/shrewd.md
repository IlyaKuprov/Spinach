# kernel/grids/shrewd.m

[MATLAB source](https://github.com/IlyaKuprov/Spinach/blob/main/kernel/grids/shrewd.m) · [Spinach Wiki](https://spindynamics.org/wiki/index.php?title=shrewd.m) · [Eden–Levitt reference](http://dx.doi.org/10.1006/jmre.1998.1427)

- Signature: `weights=shrewd(alphas,betas,gammas,max_rank,max_error)`
- `alphas`, `betas`, and `gammas`: equal-length, finite, real column vectors of active-ZYZ Euler angles in radians. The function identifies a two-angle grid when every `alpha` is exactly zero; otherwise it uses the full three-angle branch.
- `max_rank`: positive integer spherical rank. `max_error`: finite, real, non-negative scalar; the source does not impose an upper bound.
- Output `weights`: one real column entry per grid point, normalised to sum to one. These are dimensionless relative quadrature weights, not frequencies or times.

For each grid point and each rank `l=0,...,max_rank`, the function builds a matrix from Wigner D-matrix entries. With all-zero `alphas`, it keeps the single-index entries `D(l+1,l+m+1)` for `m=l,...,-l` (`(max_rank+1)^2` rows). Otherwise it keeps `D(l+m+1,l+n+1)` for all `m,n=l,...,-l` (`sum((2*l+1)^2)` rows over ranks `l=0,...,max_rank`). The right-hand side is `max_error` in every row except its first entry, which is `1-max_error`. It solves `H\v` using MATLAB backslash, takes the real part, and divides by the sum.

The resulting weights are not constrained positive by the solve: the source errors if any is exactly zero or negative, with guidance to increase `max_rank` or reduce `max_error`, respectively. The routine checks column shape, finiteness, reality, matching lengths, rank, and error scalar, but does not constrain Euler-angle ranges. For fixed inputs it has no random initialisation. It assigns orientation-grid weights; it does not calculate eigenfields, time evolution, or frequency offsets.
