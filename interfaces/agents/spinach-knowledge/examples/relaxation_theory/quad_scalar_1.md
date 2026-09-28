# examples/relaxation_theory/quad_scalar_1.m

- Signature: `quad_scalar_1()`

## Purpose

Simulates the NMR spectrum of `17O`-enriched water inside a fullerene cage. The example describes proton relaxation arising from oxygen quadrupolar relaxation together with H–O scalar coupling; the gaseous-water `17O` quadrupolar parameters are cited to [Table III](http://dx.doi.org/10.1063/1.1672122). Calculation time: seconds.

## Physical / mathematical content

The system is two protons and one `17O` nucleus at `14.1 T`. The two proton–oxygen scalar couplings are `80 Hz` each, and the proton–proton coupling is set to `0 Hz`. The oxygen has spin `5/2` and quadrupolar parameters `9.82e6` and asymmetry `0.407`. Redfield relaxation uses zero equilibrium, secular retention, and correlation time `1e-13 s`.

## Numerical / algorithmic content

The calculation uses the `sphten-liouv` basis with no approximation. It simulates a proton acquisition with `liquid` and `acquire`, using `L+` for the initial state and receiver, zero offset, a `500 Hz` sweep, `512` points, `2048` zero-fill points, and an axis in hertz. The free-induction signal is Fourier transformed with `fftshift(fft(fid,2048))` and plotted with `plot_1d`.
