# tests/kernel/test_chemical_exchange_conservation.m

Source: [tests/kernel/test_chemical_exchange_conservation.m](https://github.com/IlyaKuprov/Spinach/blob/main/tests/kernel/test_chemical_exchange_conservation.m)

## Purpose

Registered boundary regression for unsupported two-substance chemical exchange. The historical test name and manifest registration remain, but the assertion now requires `Spinach:kinetics:segmentedChemistry` and the reaction-record implementation (WP3) message rather than claiming numerical population conservation.

## Physical fixture and assertion

The fixture retains two protons at 14.1 T with zero chemical shifts, separate chemical parts `{1,2}`, equal concentrations `[1 1]`, and the symmetric exchange-rate matrix `[-3 3;3 -3]` in inverse seconds. Its `sphten-liouv` basis uses `{'none','none'}` approximations. Calling `kinetics` must raise the named boundary error; an unrelated error or silent acceptance fails the test.

There is no single-substance conservation case in this fixture. The unsupported two-site zero-column-sum assertion is deferred until segmented reaction-record chemistry is implemented. Supported single-substance flux conservation is covered separately by the kinetics generator and invariant suites.

## Inputs and outputs

```matlab
result=test_chemical_exchange_conservation()
```

No inputs. `result` is the regression result produced by `new_test_result` and updated by `test_true`.
