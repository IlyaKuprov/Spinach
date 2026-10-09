# tests/kernel/test_chemical_exchange_conservation.m

Source: [tests/kernel/test_chemical_exchange_conservation.m](https://github.com/IlyaKuprov/Spinach/blob/main/tests/kernel/test_chemical_exchange_conservation.m)

## Purpose

Registered numerical regression for two-substance chemical exchange: total concentration, total longitudinal magnetisation, matched spin transfer, and the analytic population trajectory.

## Physical fixture and assertions

Two protons at 14.1 T have zero shifts, separate chemical parts `{1,2}`, concentrations `[1 1]`, and complete spherical-tensor bases. Two directed first-order records match the protons in opposite directions at 3 s⁻¹ each. The generator conserves total unit-coordinate population and total longitudinal magnetisation; its action on the first-site longitudinal state equals three times destination minus source. A population initially `[2;0]` at 0.2 s agrees with `expm([-3 3;3 -3]*0.2)*[2;0]`. Conservation and routing use absolute and relative whole-vector tolerances of 1e-14; the trajectory uses 1e-13.

## Inputs and outputs

```matlab
result=test_chemical_exchange_conservation()
```

No inputs. `result` is the regression result produced by `new_test_result` and updated by `test_close`.
