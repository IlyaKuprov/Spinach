# kernel/cce/cce_converge.m

[MATLAB source](https://github.com/IlyaKuprov/Spinach/blob/main/kernel/cce/cce_converge.m)

## Purpose

Run independent, explicit convergence ladders and report acceptance gates.

## Behaviour

- Varies one supported field at a time from the baseline; order, bath size, graph cutoff, sample count and time-grid refinement are separate tests.
- Applies user-declared observable tolerances and work cap. Time-grid comparisons require nested points; sample-count ladders require an independent-seed check.
- Runs explicit held-out complete parameter checks when supplied; does not equate an order plateau with convergence in omitted physical dimensions.

## Inputs and outputs

```matlab
result=cce_converge(spin_system,parameters)
```

**Inputs**

- `spin_system`, `parameters` — physical model and ordinary CCE settings plus `convergence.ladders` (field and cell-row values), atol, rtol, max_work, and optional checks (cell of complete parameter structures).

**Outputs**

- `result` — baseline/run results, ladder comparisons, independent checks, work usage and convergence decision.

## References

- [Source file on GitHub](https://github.com/IlyaKuprov/Spinach/blob/main/kernel/cce/cce_converge.m)
- [cce.m](cce.md)
- [cce_diagnostics.m](cce_diagnostics.md)
