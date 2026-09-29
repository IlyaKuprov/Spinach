# kernel/cache/sle_operators.m

`[Lx,Ly,Lz,D,space_basis]=sle_operators(max_rank,int_ranks)` builds a truncated Wigner-function basis and the operators needed for lab-space rotational diffusion and interaction multiplication. `max_rank` is a positive integer. For `R = max_rank`, `space_basis` has one row `[L M N]` for every `L = 0,...,R` and every pair `M,N = L,L-1,...,-L`; its dimensions are `basis_dim`-by-3, where `basis_dim = (R+1)*(2*R+1)*(2*R+3)/3`.

`Lx`, `Ly`, and `Lz` are sparse `basis_dim`-by-`basis_dim` matrices. The raising matrix `L+` raises `M` by one at fixed `L,N`, with matrix element `sqrt(L*(L+1)-M*(M+1))`; it then sets `Lx=(L+ + L+')/2`, `Ly=(L+ - L+')/(2i)`, and `Lz=diag(M)`.

`int_ranks` is an optional row vector of distinct positive integer interaction ranks; it may be empty when only the rotation generators are needed. The validator does not impose an upper bound relating an interaction rank to `max_rank`. `D` is indexed by rank, and each requested `D{r}` is a `(2r+1)`-by-`(2r+1)` cell array. Entry `D{r}{m,n}` is a sparse `basis_dim`-by-`basis_dim` matrix for multiplication by the Wigner function with projections `M=r+1-m` and `N=r+1-n`. Its retained couplings use the product of the two Clebsch–Gordan coefficients and the normalisation factor `sqrt((2*L2+1)/(2*L+1))`. Rank 2 uses the built-in bypass formula; other ranks call `clebsch_gordan.m` through the Java virtual machine. A cached result can be loaded without that calculation.

The cache key includes `max_rank` and the requested interaction ranks in `sle_operators_rank_<max_rank>_int_<ranks>.mat`, beside the function. A cache hit loads `space_basis`, the three generators, and `D`; otherwise the results are built and a v7.3 save is attempted. A failed save warns but does not discard the returned operators. `max_rank` must be a positive integer, and `int_ranks` must be empty or a row of distinct positive integers.

[MATLAB implementation](https://github.com/IlyaKuprov/Spinach/blob/main/kernel/cache/sle_operators.m) · [Spinach Wiki](https://spindynamics.org/wiki/index.php?title=sle_operators.m)
