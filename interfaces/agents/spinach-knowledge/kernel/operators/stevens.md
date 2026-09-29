# kernel/operators/stevens.m

- Source: [kernel/operators/stevens.m](https://github.com/IlyaKuprov/Spinach/blob/main/kernel/operators/stevens.m)
- Wiki: [stevens.m on the Spin Dynamics Wiki](https://spindynamics.org/wiki/index.php?title=stevens.m)
- Signature: `S=stevens(mult,k,q)`

## Purpose

Constructs the extended Stevens operator matrix for a finite spin of multiplicity `mult`, rank `k`, and projection `q`. It uses the spin matrices returned by `pauli(mult)`; it is not a bosonic operator constructor.

## Operator construction

Set `L=pauli(mult)` and start with `S=L.p^k`. The code applies the lowering commutator `S=L.m*S-S*L.m` once for each integer from `abs(q)` through `k-1`, giving `k-abs(q)` such steps. It then symmetrises the result according to the sign of `q`:

- For `q>=0`, the result is `coeff*(S+S')/2`.
- For `q<0`, the result is `coeff*(S-S')/(2i)`.

Here `S'` is MATLAB's conjugate transpose. The input checks require a positive integer multiplicity (ultimately enforced by `pauli`), integer rank `0<=k<=12`, and integer projection `-k<=q<=k`.

## Coefficient and normalisation

The source stores coefficient vectors explicitly in `C` for ranks 0 through 12 rather than generating them from a formula. The scalar starts as `(-1)^(k-q)/C{k+1}(abs(q)+1)`; when `k` is even and `q` is odd, the code divides that scalar by an additional 2. This scalar is then used in the sign-dependent symmetrisation above. These are the implementation's normalisation conventions; no additional physical interpretation is implied here.
