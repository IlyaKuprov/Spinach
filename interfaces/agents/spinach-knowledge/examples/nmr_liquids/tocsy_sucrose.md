# examples/nmr_liquids/tocsy_sucrose.m

- Signature: `tocsy_sucrose()`

## Purpose

Simulates a liquid-state 2D TOCSY spectrum of sucrose using magnetic parameters from a vacuum DFT calculation. The source states a calculation time of seconds.

## Spin system and basis

- Parses `../standard_systems/sucrose.log` with `gparse` and converts hydrogen atoms to `1H` with `g2spinach`, using `options.min_j=1.0` and the supplied conversion parameter `31.8`.
- Sets the magnet field to `5.9`.
- Uses the `sphten-liouv` formalism with `IK-2` approximation, `scalar_couplings` connectivity, and proximity level `1`.
- Enables `greedy`, **disables `krylov`**, and sets `sys.tols.prox_cutoff=4.0`.

## Sequence and processing

- Creates the spin system and basis, then simulates with `liquid(spin_system,@tocsy,parameters,'nmr')`.
- Sets mixing time `0.100`, `lamp=1e4`, offset `800`, sweeps `[1700 1700]`, acquisition points `[512 512]`, zero filling `[2048 2048]`, observed spins `{'1H'}`, and axis units `ppm`. The initial state is proton `Lz`.
- Applies squared-cosine apodisation in both dimensions to the cosine and sine signals. The F2 transforms use the imaginary part of the cosine signal and the real part of the sine signal; these are combined as `f1_cos-1i*f1_sin` for the States signal. An F1 Fourier transform then produces the spectrum.
- Plots `abs(spectrum)` with `plot_2d`.