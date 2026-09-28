# examples/relaxation_theory/sle_esr_nitroxide_2.m

- Signature: `sle_esr_nitroxide_2()`

## Purpose

Simulates the slow-motion ESR spectrum of a nitroxide radical, as a reproduction of Figure 2 in Concilio et al. ([arXiv:1511.01667](http://arxiv.org/abs/1511.01667)). Calculation time: seconds.

## Physical / mathematical content

The system contains `14N` and an electron at `0.3343 T`. The hyperfine matrix is `1e6*gauss2mhz([10.00 3.544 11.70; 3.544 18.00 5.072; 11.70 5.072 30.00])`; the electron Zeeman matrix is `[2.0065794 -0.0007548 -0.0032848; -0.0007548 2.0056940 -0.0006008; -0.0032848 -0.0006008 2.0048920]`.

## Numerical / algorithmic content

The complete `sphten-liouv` basis is used without approximation. The SLE parameters are maximum rank `7` and correlation time `17e-9 s`. Electron `L+` is used for both initial and detection states. The simulation uses sweep `[-2.2e8,2e8]`, `1650` points and zero-fill points, `GHz-labframe` axis units, inverted axis, and first derivative; `gridfree` with `slowpass` generates the spectrum, which is plotted.
