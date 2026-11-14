# kernel/cce/cce_analyse.m

[MATLAB source](https://github.com/IlyaKuprov/Spinach/blob/main/kernel/cce/cce_analyse.m)

## Purpose

Fit an explicitly declared decay envelope without an automatic T1 claim.

## Behaviour

- Preserves input values/masks and returns status not_requested when fitting is not enabled.
- Fits A*exp(-(t/T)^beta) in a user-declared interval; exponential fixes beta=1, stretched_exponential fits it. No offset is inferred.
- Uses the selected real or magnitude envelope and explicit inverse-variance weights; optional confidence is none or linearized at the requested level. Invalid data are diagnosed rather than silently fitted.

## Inputs and outputs

```matlab
analysis=cce_analyse(result,parameters)
```

**Inputs**

- `result` — times, value and valid.
- `parameters.analysis` — fit (default false); when enabled: observable column, envelope, model, interval in seconds, positive weights, confidence and confidence_level.

**Outputs**

- `analysis` — retained data/masks, status, warnings and fit parameters/diagnostics.

## References

- [Source file on GitHub](https://github.com/IlyaKuprov/Spinach/blob/main/kernel/cce/cce_analyse.m)
- [cce_diagnostics.m](cce_diagnostics.md)
