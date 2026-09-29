# examples/relaxation_theory/quad_relaxation_1.m

[Source file](https://github.com/IlyaKuprov/Spinach/blob/main/examples/relaxation_theory/quad_relaxation_1.m) · Signature: `quad_relaxation_1()`

## Purpose and physical model

This liquid-state glycine example evaluates longitudinal and transverse relaxation rates for `14N` under quadrupolar relaxation, then compares Spinach's numerical rates with textbook analytical expressions. The source specifies a single `14N` spin and sets the magnet to 14.1. It obtains the spin quantum number from the isotope, forms the quadrupolar-interaction tensor with `eeqq2nqi` using coupling input 1.18e6 and asymmetry 0.53, and uses a correlation-time input of 1e-9. The source does not state units alongside those three parameter values.

## Calculation and output

The relaxation model is Redfield with zero equilibrium and lab-frame retention. The basis uses `sphten-liouv` with no approximation. The script builds the relaxation superoperator, evaluates the textbook rates with `rlx_nqi` using the same spin, magnet, quadrupole-coupling, asymmetry, and correlation-time inputs, and forms normalised `Lz` and `L+` states using the 2-norm.

It prints four values: longitudinal and transverse rates from the Spinach superoperator, followed by the corresponding textbook rates. The code does not include numerical output values in the source, so none are asserted here. This is a calculated comparison against analytical equations, not comparison with a measured relaxation experiment; the script prints values and does not create a plot.
