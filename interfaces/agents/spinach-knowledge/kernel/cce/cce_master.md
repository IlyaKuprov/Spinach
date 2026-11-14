# kernel/cce/cce_master.m

[MATLAB source](https://github.com/IlyaKuprov/Spinach/blob/main/kernel/cce/cce_master.m)

## Purpose

Propagate open-system conditional or generalized cluster dynamics.

## Behaviour

- Optional public fields are checked by `cce_validate`; unknown, reserved, and method-inapplicable requests are rejected rather than silently ignored.

- Selective explicit-generator pulses use `cce_hamiltonian` metadata dimensions to embed the target in the retained central-first tensor order.

- Uses the local GKSL generator and the compiled ideal-pulse schedule for `mecce` or `megcce`.
- Conditional blocks produce branch coherence; generalized density propagation supplies the same raw contract as the unitary evaluator.
- Checks local trace, Hermiticity and positivity diagnostics where defined; product-truncated populations are not a universal relaxation-time method.

- Generalized coherence is raw by default, matching `gcce`; only an explicit `normalized=true` divides by the initial central coherence.

- Explicit operator pulses target the complete ordered central register or one physical site; unsupported multi-site targets are rejected instead of ignored.

- Stock propagation uses a local copy with `prop_chop=0`: generator, Taylor-term and squaring clean-up must not discard dynamics outside the reported numerical bound. The caller's tolerance policy is unchanged.

## Inputs and outputs

```matlab
raw=cce_master(model,cluster,context,schedule,parameters)
```

**Inputs**

- `model`, `cluster`, `context`, `schedule`, `parameters` — local physical/state inputs, method, ideal sequence, explicit dissipation jumps and numerical limits.

**Outputs**

- `raw` — cluster observable/density, validity and numerical error metadata for `cce_reduce`.

## References

- [Source file on GitHub](https://github.com/IlyaKuprov/Spinach/blob/main/kernel/cce/cce_master.m)
- [cce_liouvillian.m](cce_liouvillian.md)
- [cce_observe.m](cce_observe.md)
- [Method reference](https://doi.org/10.1103/PhysRevLett.132.250401)
