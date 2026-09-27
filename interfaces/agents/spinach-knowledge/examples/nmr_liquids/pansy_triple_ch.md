# examples/nmr_liquids/pansy_triple_ch.m

- Signature: `pansy_triple_ch()`

## Purpose

Triple-channel PANSY experiment on glycine with natural content of 13C isotope. Coordinates, shieldings, and J-couplings were computed with DFT. Calculation time: seconds.

## Implementation

- Read the vacuum-DFT glycine system from `../standard_systems/glycine.log` using `gparse` and `g2spinach`. Map H, C, and N to `1H`, `13C`, and `15N` with values `[31.5 189.2 400.3]`; set `options.min_j=3.0`, `options.no_xyz=1`, and `sys.magnet=5.9`.
- Create the spin system with a `sphten-liouv` basis, `IK-2` approximation, `scalar_couplings` connectivity, and proximity level `1`.
- Set channel sweeps to `[1700 4500 8000]`, offsets to `[900 3200 4000]`, acquisition points to `[128 128 128]`, zero filling to `[512 512 512]`, spins to `{'1H','13C','15N'}`, and axis units to `ppm`.
- Generate isotopomers by diluting first on `13C`, then on `15N`. For each isotopomer, build the basis and run `liquid(subsystem,@pansy_triple,parameters,'nmr')`.
- Apply two-dimensional `sqcos` apodisation to `fid.aa`, `fid.ab`, and `fid.ac`; zero-fill, Fourier-transform with `fft2`, shift with `fftshift`, and sum each result across isotopomers.
- Plot the real spectra as three 2D blocks: `1H`–`1H`, `1H`–`13C`, and `1H`–`15N`, using their corresponding saved channel axes.