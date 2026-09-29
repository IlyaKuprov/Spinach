# examples/nmr_liquids/coloc_test.m

- MATLAB implementation: [examples/nmr_liquids/coloc_test.m](https://github.com/IlyaKuprov/Spinach/blob/main/examples/nmr_liquids/coloc_test.m)

- Signature: `coloc_test()`

## Purpose

A compact two-spin 1H-13C COLOC simulation that illustrates long-range scalar-coupling transfer. The driver describes a calculation time of seconds.

## Spin system and transfer

The two spins are 1H and 13C at 11.7 T, with isotropic shifts 4.0 and 75.0 ppm and a scalar coupling of 5.0 Hz. The source also sets the carbon diagonal coupling entry to zero. It builds the full sphten-Liouville basis (approximation 'none') and calls the COLOC pulse program in NMR mode.

In the pulse program the initial state is proton longitudinal magnetisation. A proton 90-degree pulse begins the sequence; the indirect evolution is embedded in an echo with simultaneous 180-degree proton and carbon pulses. The program selects -1 proton coherence, applies proton and carbon 90-degree pulses, evolves for delta2, selects +1 carbon coherence, decouples the proton, and acquires on a carbon L+ detection operator. The driver sets delta2 = 30 ms. It leaves delta1 unset; the pulse program defaults delta1 to half the maximum F1 evolution time implied by the sweep (25.5 ms for these F1 settings). The cited implementation follows Fig. 1b of the COLOC sequence and omits the dashed pulses during delta2.

## Acquisition and display

The parameter order is F1 = 1H and F2 = 13C. Offsets are [2250, 5000] Hz, sweep widths [5000, 12000] Hz, and the acquired grid is [256, 256] points, zero-filled to [512, 512]. The axes are labelled in ppm. Cosine apodisation is applied in both dimensions, followed by a shifted two-dimensional Fourier transform; the plot shows the magnitude spectrum.

The pulse-sequence implementation cites the original COLOC paper: [DOI 10.1016/0022-2364(84)90136-7](https://doi.org/10.1016/0022-2364(84)90136-7).

## Scope

The example contains only the stated two-spin system and a simulated pulse sequence; it does not model a larger molecule or report a comparison with measured data. The magnitude display suppresses the sign and phase information of the complex transformed signal.
