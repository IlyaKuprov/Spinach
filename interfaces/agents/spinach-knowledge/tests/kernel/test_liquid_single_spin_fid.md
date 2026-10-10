# tests/kernel/test_liquid_single_spin_fid.m

**Source:** [test_liquid_single_spin_fid.m](https://github.com/IlyaKuprov/Spinach/blob/main/tests/kernel/test_liquid_single_spin_fid.m)

## Purpose

Checks that a zero-offset isolated proton with no relaxation has constant transverse magnetisation throughout a liquid-state free induction decay. The production liquid/acquire pathway must preserve its initial detected coherence.

## Numerical specification

The complete one-proton spherical-tensor Liouville basis is used, with zero chemical shift and no spin-spin coupling. The eight-point acquisition has a 1000 Hz sweep. Every FID point is compared with the first point, using absolute and relative tolerances of `1e-12`.

The proton carrier angular frequency is fixed at `3772062842.904` rad/s, the original 14.1 T regression value. Dividing that frequency by `spin('1H')` determines the field, so a literature magnetic-moment update does not change the numerical model. This test checks signal conservation, not the adopted value of the proton moment.

## Inputs and outputs

```matlab
result=test_liquid_single_spin_fid()
```

No inputs. The output is a regression result structure containing the constant-FID comparison and its explanatory messages.
