# examples/relaxation_theory/sle_esr_nitroxide_2.m

[MATLAB source](https://github.com/IlyaKuprov/Spinach/blob/main/examples/relaxation_theory/sle_esr_nitroxide_2.m)

## Purpose

This example calculates a slow-motion ESR spectrum for a nitroxide radical. The source describes the calculation as a reproduction of Figure 2 in Concilio et al. ([arXiv:1511.01667](https://arxiv.org/abs/1511.01667)); the plotted trace is a simulation, not an experimental spectrum supplied by the example. The source lists a calculation time of seconds.

## Spin model and motion

The spin system is one 14N nucleus and one electron; the source sets `sys.magnet=0.3343`. Its 14N-electron hyperfine input is the matrix

```
[10.00  3.544 11.70; 3.544 18.00 5.072; 11.70 5.072 30.00]
```

passed through `gauss2mhz`. The source does not state a separate unit for the tensor entries. The electron Zeeman tensor is set to `[2.0065794 -0.0007548 -0.0032848; -0.0007548 2.0056940 -0.0006008; -0.0032848 -0.0006008 2.0048920]`; no unit is specified for this matrix in the file. The basis is `sphten-liouv` with no basis approximation.

The Stochastic Liouville Equation (SLE) calculation uses maximum rank 7 and `parameters.tau_c=17e-9`. Both the initial state and detection coil are the electron raising state, `L+`; the selected spin is `E`, and the decoupling list is empty. The spectrum is generated with `gridfree` and `slowpass` in ESR mode; this script does not set a separate Bloch-Redfield relaxation model.

## Spectrum and plot

The sweep is `[-2.2e8, 2e8]` on the `GHz-labframe` axis, with 1650 points and 1650 zero-fill points. The axis is inverted and the first derivative is requested. The script plots the real part of the calculated signal with `plot_1d`; the example specifies no measured line positions or intensities.
