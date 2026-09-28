# examples/dnp_sol/cross_effect_field_scan_1.m

- Signature: `cross_effect_field_scan_1()`

## Purpose

Calculates the steady-state proton magnetisation under microwave irradiation for a powder-averaged cross-effect DNP system as a function of magnetic-field offset. The source estimates hours to run.

## Model

At 18.78 T, the four-spin system is `E,E,14N,1H`. It includes anisotropic g tensors for both electrons, a `^14N` quadrupolar tensor, electron–nitrogen and electron–electron couplings, and the specified spin coordinates. Relaxation uses the `t1_t2` model with diagonal retention, zero equilibrium, and temperature 10 K. The full sphten-Liouville basis is used.

## Calculation

The script sets 10 MHz microwave power, sweeps 256 field offsets from −0.08 to +0.04 T, and powder-averages on `rep_2ang_1600pts_sph`. It calls `powder` with `dnp_field_scan` in ESR mode and plots the real proton `Lz` expectation against magnetic-field offset. The configured powder method is `backslash`.
