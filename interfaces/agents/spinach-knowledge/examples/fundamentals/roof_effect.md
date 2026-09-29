# examples/fundamentals/roof_effect.m

- MATLAB implementation: [examples/fundamentals/roof_effect.m](https://github.com/IlyaKuprov/Spinach/blob/main/examples/fundamentals/roof_effect.m)

## Purpose

Illustrate the roof effect in the spectrum of a strongly J-coupled two-proton system while bringing the two resonance offsets closer together.

## System and acquisition

The source configures two `1H` spins at 5.9 T, scalar Zeeman values 0.95 and 1.45, and a 7.0 Hz scalar coupling. It uses the `sphten-liouv` basis without approximation. For each plotted case, the source replaces the two Zeeman matrices with values corresponding to offsets around 1.2 ppm: `1.2 + ppm` and `1.2 - ppm`, where `ppm` takes 0.2, 0.05, 0.0125, and 0.00625. Thus the sweep parameter controls the symmetric separation; the initially configured scalar shifts are not used unchanged in those four simulations.

The acquisition uses proton `L+` for both the initial state and receiver, runs `liquid` with the NMR assumption, and sets offset and sweep to 300 Hz with 1,024 points. The FID is exponentially apodised with parameter 10, zero-filled to 4,096 points, Fourier transformed, and shifted for plotting. The plotted ordinate is the real spectrum and the axis is in Hz with inversion enabled.

## Output and limitations

The output is a four-panel set of spectra for the four offset settings. The source specifies no numerical roof-effect metric, acceptance threshold, or automated assertion; the page therefore does not claim that a particular spectral shape or intensity ratio was measured. The source contains no cited external reference for this example.

[Source example](../../../../../examples/fundamentals/roof_effect.m).
