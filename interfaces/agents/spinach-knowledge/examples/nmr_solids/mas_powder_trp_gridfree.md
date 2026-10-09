# examples/nmr_solids/mas_powder_trp_gridfree.m

- MATLAB implementation: [examples/nmr_solids/mas_powder_trp_gridfree.m](https://github.com/IlyaKuprov/Spinach/blob/main/examples/nmr_solids/mas_powder_trp_gridfree.m)

Source: [examples/nmr_solids/mas_powder_trp_gridfree.m](../../../../../examples/nmr_solids/mas_powder_trp_gridfree.m)

- Signature: `mas_powder_trp_gridfree()`

## Model and data

This example calculates the `13C` MAS spectrum of tryptophan powder assuming `1H` decoupling. It imports tryptophan data from `trp_xray.out`, mapping C and N to `13C` and `15N`, and sets 9.4 T. The source describes the coordinates as X-ray data, chemical-shift anisotropies as DFT estimates, and isotropic shifts as experimental inputs. These assigned shifts are model inputs, not calculated peak positions.

For the two unit-cell molecules, the first eight indexed shift values are 124.2, 110.1, 114.7, 118.0, 119.3, 107.5, 134.9, and 125.0 ppm. The final three are 26.8, 54.6, and 174.4 ppm for the first molecule, and 28.0, 52.1, and 173.3 ppm for the second. The source cites [the paper at DOI 10.1126/sciadv.aaw8962](https://doi.org/10.1126/sciadv.aaw8962) for further particulars on the polyadic representation.

## Grid-free Fokker-Planck calculation

The basis uses `sphten-liouv`, `IK-0`, longitudinal `15N`, projection +1, and interaction level 3. The rotor rate is 14 kHz about axis vector `[1 1 1]`; rank 11 is used without a named orientation-grid parameter. The first system enables `greedy` and `polyadic`; the second adds `gpu`. Each conformation is propagated by `gridfree` in NMR mode with the `nmr` assumption, and its FID is added to the first. Acquisition is on `13C` with `L+` initial and receiver states, 100 kHz sweep, 2048 points, 8192-point zero filling, and zero offset; the display axis is ppm and inverted.

## Spectrum processing

The summed FID is exponentially apodised with parameter 6 and Fourier transformed. The plotted trace is the real part of the calculated spectrum.
