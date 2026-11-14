# kernel/cce/cce_noise.m

[MATLAB source](https://github.com/IlyaKuprov/Spinach/blob/main/kernel/cce/cce_noise.m)

## Purpose

Evaluate bath-noise autocorrelation for additive connected-cluster expansion.

## Behaviour

- Builds the sum of noise terms completely supported inside the cluster; refuses a split multi-site term.
- Propagates with the bath Hamiltonian after tracing out the central register at infinite temperature and applies the requested ordered or symmetrized correlation convention.
- Optional connected subtraction removes the mean-product contribution; the result retains units squared of the supplied observable operators.

## Inputs and outputs

```matlab
raw=cce_noise(model,cluster,context,schedule,parameters)
```

**Inputs**

- `model`, `cluster`, `context`, `schedule`, `parameters` — cluster inputs; `noise.ordering=ordered|symmetrized`, logical connected, and terms with ordered bath support and operator matrices.

**Outputs**

- `raw` — time-dependent autocorrelation with the standard cluster IDs, validity and error fields; additive reduction is selected by method noise.

## References

- [Source file on GitHub](https://github.com/IlyaKuprov/Spinach/blob/main/kernel/cce/cce_noise.m)
- [cce_reduce.m](cce_reduce.md)
- [cce_spectrum.m](cce_spectrum.md)
