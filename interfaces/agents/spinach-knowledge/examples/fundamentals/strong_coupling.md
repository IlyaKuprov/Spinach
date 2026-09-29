# examples/fundamentals/strong_coupling.m

- Signature: `strong_coupling()`
- Source: [`examples/fundamentals/strong_coupling.m`](https://github.com/IlyaKuprov/Spinach/blob/main/examples/fundamentals/strong_coupling.m)

## Purpose

A compact liquid-state NMR demonstration of a strongly coupled two-proton system, from spin-system setup through acquisition and spectrum plotting.

## Spin model and basis

The source sets two `1H` spins, a magnet-field parameter of `5.9`, Zeeman scalar parameters `1.0` and `1.5`, and a single scalar coupling of `7.0`. It uses the full `sphten-liouv` basis with `approximation='none'`; the source does not spell out units for the field, Zeeman entries, or coupling.

## Acquisition and processing

The initial state and receiver are both the `1H` `L+` state. The acquisition uses an empty decoupling list, offset `300`, sweep width `300`, `1024` points, and zero-filling to `4096`; the plotted axis is explicitly in `Hz` and is inverted. The source calls `liquid` with the `acquire` callback in NMR mode, applies exponential apodisation with parameter `10`, Fourier-transforms the FID, and plots the real part of the shifted spectrum.

## Output and limits

The result is a plotted, processed 1D spectrum; the source does not provide an expected numerical spectrum or assertion threshold. The acquisition sweep and displayed axis use hertz as specified by `axis_units='Hz'`; other units are not stated in the file.
