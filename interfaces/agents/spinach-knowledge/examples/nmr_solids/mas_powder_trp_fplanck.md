# examples/nmr_solids/mas_powder_trp_fplanck.m

- MATLAB implementation: [examples/nmr_solids/mas_powder_trp_fplanck.m](https://github.com/IlyaKuprov/Spinach/blob/main/examples/nmr_solids/mas_powder_trp_fplanck.m)

Source: [examples/nmr_solids/mas_powder_trp_fplanck.m](../../../../../examples/nmr_solids/mas_powder_trp_fplanck.m)

- Signature: `mas_powder_trp_fplanck()`

## Model and data

This example calculates the `13C` MAS spectrum of tryptophan powder assuming `1H` decoupling, using the Fokker-Planck MAS formalism. It imports tryptophan data from `trp_xray.out`, maps C and N to `13C` and `15N`, and sets 9.4 T. The source identifies X-ray coordinates and DFT-estimated chemical-shift anisotropies; its indexed isotropic shifts are described as experimental inputs, not simulated output.

For the two unit-cell molecules, the first eight indexed shift values are 124.2, 110.1, 114.7, 118.0, 119.3, 107.5, 134.9, and 125.0 ppm. The first molecule uses 26.8, 54.6, and 174.4 ppm for the final three; the second uses 28.0, 52.1, and 173.3 ppm. These are Zeeman-tensor inputs, not reported spectrum positions.

## Fokker-Planck calculation and acquisition

The basis is `sphten-liouv` with `IK-0`, longitudinal `15N`, projection +1, and interaction level 3. The rotor rate is 14 kHz about axis vector `[1 1 1]`; the calculation uses rank 11 and the named `leb_2ang_rank_11` grid. Each system is propagated by `singlerot` with `acquire` in NMR mode, and the two FIDs are added. Acquisition uses a 100 kHz sweep, 2048 points, 8192-point zero filling, and `13C` `L+` initial and receiver states. The source enables the GPU for the second molecule; the GPU-enable line for the first is commented out.

## Spectrum processing

The summed FID is exponentially apodised with parameter 6, Fourier transformed, and plotted as its real spectrum. The simulation settings are not measured spectral output.
