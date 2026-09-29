# examples/relaxation_theory/quad_scalar_1.m

[Source file](https://github.com/IlyaKuprov/Spinach/blob/main/examples/relaxation_theory/quad_scalar_1.m) · Signature: `quad_scalar_1()`

## Purpose

Calculates an NMR spectrum for `17O`-enriched water inside a fullerene cage. The source describes proton relaxation from the combination of oxygen quadrupolar relaxation and H–O scalar coupling; it cites the quadrupolar parameters for gaseous (assumed) water in [Table III of the referenced paper](https://doi.org/10.1063/1.1672122). This is a calculated spectrum, not a measured spectrum. The example comments give a calculation time of seconds.

## Model and parameters

The spin list is two `1H` spins and one `17O` spin, with `sys.magnet=14.1`. The proton–proton scalar entry is zero; each proton–oxygen scalar entry is `80`. The oxygen quadrupolar matrix is built by `eeqq2nqi(9.82e6,0.407,5/2,[0 0 0])`. These numerical interaction inputs are reproduced as written; the source does not label their units. Relaxation is set to `redfield`, with zero equilibrium, secular retention, and `tau_c=1e-13` (unit not stated). The basis is the full `sphten-liouv` basis with no approximation.

## Acquisition and output

The acquisition uses `1H` for both `rho0` and the coil, no decoupling, zero offset, sweep `500`, `512` points, zero-fill `2048`, and axis units `Hz`. `liquid(...,'nmr')` produces the FID; the script applies `fftshift(fft(fid,parameters.zerofill))` and plots its real part with `plot_1d`. The source specifies the initial and detection operators as `L+` on `1H`.
