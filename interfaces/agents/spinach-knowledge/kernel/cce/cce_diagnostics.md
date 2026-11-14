# kernel/cce/cce_diagnostics.m

[MATLAB source](https://github.com/IlyaKuprov/Spinach/blob/main/kernel/cce/cce_diagnostics.m)

## Purpose

Report validity, conditioning, sampling and central-density diagnostics.

## Behaviour

- Constructs per-observable invalid reasons, valid-time masks/fraction and numerical budget checks from result values, masks and error bounds.
- Reports cluster cost/order information and sample statistics when available; density outputs are checked for trace, Hermiticity and minimum eigenvalue.
- Never clips observables or projects densities to hide failed validity or physicality checks.

## Inputs and outputs

```matlab
diagnostics=cce_diagnostics(model,raw,result,parameters)
```

**Inputs**

- `result` — CCE value, valid mask, abs_error_bound and retained metadata.
- `parameters` — numerical tolerances and observable declarations.

**Outputs**

- `diagnostics` — numerical/conditioning, cluster, sample and density diagnostic fields.

## References

- [Source file on GitHub](https://github.com/IlyaKuprov/Spinach/blob/main/kernel/cce/cce_diagnostics.m)
- [cce.m](cce.md)
- [cce_reduce.m](cce_reduce.md)
