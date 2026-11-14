# kernel/cce/cce_cluster.m

[MATLAB source](https://github.com/IlyaKuprov/Spinach/blob/main/kernel/cce/cce_cluster.m)

## Purpose

Propagate one unitary conditional or generalized cluster.

## Behaviour

- Conditional CCE evaluates branch traces `Tr(Ua*rhoB*Ub^dagger)`; generalized CCE propagates the complete central-register-plus-cluster density.
- Uses the same exterior context for every cluster and returns complex values with numerical validity/error information.
- Generalized central densities are passed through `cce_observe`; static Hermitian eigensystems supply exact local exponentials and measured backward residuals.

## Inputs and outputs

```matlab
raw=cce_cluster(model,cluster,context,schedule,parameters)
```

**Inputs**

- `model`, `cluster`, `context`, `schedule`, `parameters` — prepared model; record with `id` and `sites` (`id=0`, empty sites for baseline); one state context; compiled sequence; method, observable and numerics.

**Outputs**

- `raw` — cluster/state/configuration IDs, `value`, `valid`, `abs_error_bound`, central density where applicable, and propagation diagnostics.

## References

- [Source file on GitHub](https://github.com/IlyaKuprov/Spinach/blob/main/kernel/cce/cce_cluster.m)
- [cce_hamiltonian.m](cce_hamiltonian.md)
- [cce_observe.m](cce_observe.md)
- [cce_reduce.m](cce_reduce.md)
