# kernel/multiprop.m

- Signature: `rho=multiprop(spin_system,P,rho,N)`

## Meaning

Applies a square propagator `P` exactly `N` times without forming `P^N` as a separate matrix. Binary exponentiation reduces the number of matrix products; a validated `N=0` returns the original `rho` unchanged. Squaring uses panel multiplication with `spin_system.tols.prop_chop`; ordinary cleanup controls storage, and numerical agreement within roundoff—not bitwise identity—is the contract.

For `sphten-liouv`, `zeeman-liouv`, and `zeeman-wavef`, each active power acts as `rho=P*rho`. In `zeeman-hilb`, `rho=P*rho*P'`, where MATLAB `P'` is the conjugate transpose. The output retains `rho`'s dimensions.

## Inputs and guards

- `spin_system.bas.formalism` must be one of the four formalisms above, and `spin_system.tols.prop_chop` must be a non-negative real scalar.
- `P` must be a finite numeric square matrix; `rho` a finite numeric matrix.
- `N` must be a real numeric scalar representing a non-negative integer. Non-integer MATLAB numeric classes are also checked for finiteness and an upper bound of `flintmax`; integer classes are checked for non-negativity. The routine converts accepted values to `uint64` for bit processing.
- Dimensions must agree. For the Hilbert-space density-matrix formalism, `rho` must be square; in the other formalisms its row count must match `P`.

## References

- [MATLAB source](https://github.com/IlyaKuprov/Spinach/blob/main/kernel/multiprop.m)
- [Spin Dynamics Wiki](https://spindynamics.org/wiki/index.php?title=multiprop.m)
