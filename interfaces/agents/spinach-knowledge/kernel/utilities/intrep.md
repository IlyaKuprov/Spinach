# kernel/utilities/intrep.m

**Source:** [https://github.com/IlyaKuprov/Spinach/blob/main/kernel/utilities/intrep.m](https://github.com/IlyaKuprov/Spinach/blob/main/kernel/utilities/intrep.m)

## Purpose

Transforms a laboratory-frame Hamiltonian into the interaction representation with respect to a specified Hamiltonian, to a specified order in perturbation theory, following the auxiliary matrix method described in [https://doi.org/10.1063/1.4928978](https://doi.org/10.1063/1.4928978).

## Behaviour

- Syntax: `Hr=intrep(spin_system,H0,H,T,order)`.
- Validates consistency: `H` and `H0` must be Hermitian; `order` must be a non-negative integer or `Inf`.
- Confirms that `T` is a period of the `H0` propagator by computing `P=propagator(spin_system,H0,T)` and checking `norm(P-speye(size(P)),1)>1e-6`; an error is raised if the check fails.
- Computes and reports the 1-norms of `H0` (`norm(H0,1)`) and `H1=H-H0` (`norm(H-H0,1)`), and reports the rotating frame period in seconds.
- Errors if `norm_h1>norm_h0`, since `H1` must be a perturbation.
- Order handling:
  - `order=0`: shortcut for high field, `Hr=H-H0`.
  - `order=Inf`: shortcut for infinite order, `Hr=(1i/T)*logm(expm(full(-1i*H*T)))`.
  - Otherwise: computes derivatives via `dirdiff(spin_system,H0,H-H0,T,order+1)`, forms the first term `Hr=(1i/T)*(D{1}'*D{2})`, and adds the remaining series terms `Hr=Hr+(1i/T)*nchoosek(n-1,k-1)*D{n-k+1}'*D{k+1}/factorial(n)` for `n=2:order`, `k=1:n`.
- Symmetrises and cleans the output: `Hr=clean_up(spin_system,(Hr+Hr')/2,spin_system.tols.liouv_zero)`.
- Reports matrix density statistics: dimension, number of nonzeros (`nnz`), density percentage, and sparsity flag (`issparse`).
- The source notes that the auxiliary matrix method is massively faster than either commutator series or diagonalisation.

## Inputs and outputs

**Inputs**

- `spin_system` — Spinach spin system object.
- `H0` — Hamiltonian with respect to which the interaction representation transformation is done, typically the Zeeman Hamiltonian.
- `H` — laboratory frame Hamiltonian `H0+H1` to be transformed, typically the full Hamiltonian.
- `T` — period of the `H0` propagator.
- `order` — perturbation theory order in the rotating frame transformation; may be `Inf`.

**Outputs**

- `Hr` — Hamiltonian in the interaction representation.

## References

- [https://doi.org/10.1063/1.4928978](https://doi.org/10.1063/1.4928978)
- [Spinach Wiki: intrep.m](https://spindynamics.org/wiki/index.php?title=intrep.m)
