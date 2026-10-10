# kernel/conventions/transforms/ham2nqi.m

## Purpose

Decomposes a single-spin Hamiltonian in the Zeeman basis into its Zeeman vector and quadrupolar coupling tensor.

## Signature

`[omega,Q]=ham2nqi(H)`

## Decomposition convention

The returned quantities reconstruct the Hamiltonian as:

`H = omega(1)*Sx + omega(2)*Sy + omega(3)*Sz + [Sx Sy Sz]*Q*[Sx Sy Sz].'`

Here `Sx`, `Sy`, and `Sz` are the Cartesian operators from `pauli(mult)`; the tensor components are expressed in that same Cartesian operator basis. The function does not rotate the basis. Both `omega` and `Q` are in rad/s.

The code extracts each Zeeman component as `trace(Si'*H)/norm(Si,'fro')^2` for `Si = S.x, S.y, S.z`, then takes the real part. For multiplicity greater than two it obtains the rank-2 coefficients from `T{5}` through `T{9}` and converts them with `sphten2mat([],[],rank2)`; for a 2x2 spin-1/2 Hamiltonian it returns `Q=zeros(3,3)`.

## Input and outputs

- `H`: numeric square matrix in the Zeeman basis. It must be Hermitian, have `abs(trace(H)) <= eps()*norm(H,2)`, and contain at least four elements (at least 2x2). The implementation also reconstructs the Hamiltonian and errors if `norm(H-HR,2) > 1e-6*norm(H,2)`, identifying terms beyond the supported linear and quadratic form.
- `omega`: 1x3 real vector of Larmor-frequency components in rad/s.
- `Q`: 3x3 real symmetric traceless quadrupolar coupling tensor in rad/s.

## References

- MATLAB source: [kernel/conventions/transforms/ham2nqi.m](https://github.com/IlyaKuprov/Spinach/blob/main/kernel/conventions/transforms/ham2nqi.m)
- [Spinach Wiki](https://spindynamics.org/wiki/index.php?title=ham2nqi.m)
