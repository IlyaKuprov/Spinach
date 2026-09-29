# examples/kinetics/exchange_asymmetric.m

## Purpose and callable context

The no-argument MATLAB function is a compact two-spin asymmetric chemical-exchange example. The source describes a calculation time of seconds. Its numeric settings are given below as written; the source does not annotate units for the magnet field, scalar offsets, kinetic rates, acquisition offset/sweep, or exponential apodisation parameter.

## Physical model and parameters

Two 1H spin environments are assigned scalar Zeeman values 0 and 3, with exchange parts `{1,2}`. The rate matrix is `[-500 2000; 500 -2000]`, and the specified concentration weights are `[2000 500]`; the unequal off-diagonal rates make the exchange asymmetric. The source sets `sys.magnet=14.1`. It constructs the `sphten-liouv` basis with approximation `none`.

## Acquisition and observable

The initial density operator is the chemical-state-aware 1H raising operator, and the coil operator is the ordinary 1H raising operator. With no decoupling, acquisition calls `liquid` using `@acquire` in NMR mode. The parameters are offset 900, sweep 5000, 512 acquired points and 1024 zero-filled points; the frequency axis is labelled ppm and inverted. The FID is exponentially apodised with parameter 6, Fourier-transformed, and the real part of the shifted spectrum is plotted with `plot_1d`. The source contains no numeric spectrum or measured peak positions, so none are asserted here.

[MATLAB source](https://github.com/IlyaKuprov/Spinach/blob/main/examples/kinetics/exchange_asymmetric.m)
