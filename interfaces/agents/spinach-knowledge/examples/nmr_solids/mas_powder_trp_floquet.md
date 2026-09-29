# examples/nmr_solids/mas_powder_trp_floquet.m

- MATLAB implementation: [examples/nmr_solids/mas_powder_trp_floquet.m](https://github.com/IlyaKuprov/Spinach/blob/main/examples/nmr_solids/mas_powder_trp_floquet.m)

Source: [examples/nmr_solids/mas_powder_trp_floquet.m](../../../../../examples/nmr_solids/mas_powder_trp_floquet.m)

- Signature: `mas_powder_trp_floquet()`

## Model and data

This example calculates the `13C` MAS spectrum of tryptophan powder under the stated assumption of `1H` decoupling. It imports tryptophan data from `trp_xray.out`, maps C and N to `13C` and `15N`, and sets the field to 9.4 T. The source describes the coordinates as X-ray data, chemical-shift anisotropies as DFT estimates, and isotropic shifts as experimental inputs. The spectrum is calculated from those inputs; it is not a reported measured spectrum.

The code separately builds the two unit-cell molecules. Its 11 indexed isotropic-shift assignments are 124.2, 110.1, 114.7, 118.0, 119.3, 107.5, 134.9, 125.0, 26.8, 54.6, and 174.4 ppm for the first; the second retains the first eight values and uses 28.0, 52.1, and 173.3 ppm for the last three. These Zeeman-tensor inputs are not calculated peak positions.

## Floquet calculation and acquisition

Both systems use the `sphten-liouv` formalism with the `IK-0` approximation, longitudinal `15N`, projection +1, and interaction level 3. The setup uses a 14 kHz rotor rate with axis vector `[1 1 1]`, rank 11, and the `leb_2ang_rank_11` orientation grid. For each molecule, the script calls `floquet` in NMR mode and adds its FID to the other. Acquisition is on `13C` with `L+` initial and receiver states, a 100 kHz sweep, 2048 points, 8192-point zero filling, and zero offset. The source sets the display axis to ppm and inverts it.

## Spectrum processing

The summed FID receives exponential apodisation with parameter 6, then a Fourier transform. The plotted trace is the real part of the calculated spectrum.
