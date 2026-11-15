# tests/kernel/test_cce_examples.m

[MATLAB source](https://github.com/IlyaKuprov/Spinach/blob/main/tests/kernel/test_cce_examples.m)

## Purpose

Run all nineteen CCE application examples and check their independent small-system references.

## Behaviour

- Executes the CCE, gCCE, ME-CCE and partitioned-CCE examples with figures hidden, then restores the original figure-visibility setting.
- Checks finite numerical results, validity and independently evaluated reference errors returned by the examples.
- Uses observed-spin density matrices and selected-transition signals with the same contracts as the CCE kernel.
- The reduced application models illustrate algorithms; they are not material fits or converged powder calculations.

## Calling convention

```matlab
result=test_cce_examples()
```

**Inputs**

- None; run with Spinach and its test libraries on the MATLAB path.

**Outputs**

- `result` — standard test-harness record containing checks, failure messages and diagnostics.

## References

- [Registered tests](https://github.com/IlyaKuprov/Spinach/blob/main/tests/lib/test_manifest.m)
- [CCE application examples](../../examples/cce/)

The DEER example check also excludes the pump resonance and requires an empty selected event, guarding frequency selection rather than a fixed physical target.
