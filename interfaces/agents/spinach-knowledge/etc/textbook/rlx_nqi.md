# etc/textbook/rlx_nqi.m

## Signature

`[r1,r2,t1,t2]=rlx_nqi(I,omega,C_q,eta_q,tau_c)`

## Purpose

Calculates Redfield longitudinal and transverse relaxation rates for a nucleus with spin `I` subject to quadrupolar relaxation in an isotropically tumbling liquid. The implementation accepts integer or half-integer `Ige1`.

## Model and calculation

The routine constructs the quadrupolar tensor from the coupling constant and asymmetry, evaluates its rank-2 Blicharski invariant, and applies the spectral-density terms at zero, the nuclear Zeeman frequency, and twice that frequency. It then returns relaxation times as the reciprocals of the corresponding rates.

## Inputs

- `I`: nuclear spin quantum number (integer or half-integer, at least 1).
- `omega`: nuclear Zeeman frequency in radians per second.
- `C_q`: quadrupolar coupling constant `e^2 q Q/h` in Hz.
- `eta_q`: quadrupolar tensor asymmetry.
- `tau_c`: positive rotational correlation time in seconds.

## Outputs

- `r1`: longitudinal relaxation rate in Hz.
- `r2`: transverse relaxation rate in Hz.
- `t1`: longitudinal relaxation time in seconds, `1/r1`.
- `t2`: transverse relaxation time in seconds, `1/r2`.
