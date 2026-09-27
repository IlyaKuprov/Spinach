# examples/relaxation_theory/sle_solid_limit.m

- Signature: `sle_solid_limit()`

## Purpose

Calculates nitroxide ESR spectra with the stochastic Liouville equation across increasing rotational correlation times to illustrate the solid-limit regime. Calculation time: hours.

## Physical / mathematical content

The system is `14N` plus an electron at `0.3343 T`. The hyperfine tensor is `1e6*gauss2mhz([10.00 3.544 11.70; 3.544 18.00 5.072; 11.70 5.072 30.00])`; the electron Zeeman tensor is `[2.0065794 -0.0007548 -0.0032848; -0.0007548 2.0056940 -0.0006008; -0.0032848 -0.0006008 2.0048920]`.

## Numerical / algorithmic content

The complete `sphten-liouv` basis is used without approximation. Electron `L+` supplies the initial and detection states. The sweep is `[-2.2e8,2e8]` with `240` points and zero-fill points, `GHz-labframe` units, inverted axis, and no derivative. Four SLE spectra are generated with `gridfree` and `slowpass` for rank/correlation-time pairs `(3,1e-9 s)`, `(7,1e-8 s)`, `(15,1e-7 s)`, and `(30,1e-6 s)`, and plotted in a row.
