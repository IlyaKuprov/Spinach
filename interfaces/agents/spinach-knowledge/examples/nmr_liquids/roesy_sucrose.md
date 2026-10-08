# examples/nmr_liquids/roesy_sucrose.m

[Spinach source](https://github.com/IlyaKuprov/Spinach/blob/main/examples/nmr_liquids/roesy_sucrose.m)

This wrapper simulates a proton ROESY spectrum for sucrose using molecular magnetic parameters read from a vacuum-DFT log, not an experimental spectrum. It parses `../standard_systems/sucrose.log` with `gparse` and `g2spinach`, selecting `1H` and setting `options.min_j=1.0`. The source header says the magnetic parameters were computed with DFT; the source comment estimates minutes of calculation time, but no calculation or runtime measurement was performed for this draft.

At `sys.magnet=5.9` (5.9 T by Spinach convention), the spin system uses a `sphten-liouv` basis, `IK-2` approximation, scalar-coupling connectivity, and proximity level 3. The wrapper specifies Redfield relaxation, zero equilibrium, secular relaxation terms, and correlation time `200e-12` s (200 ps). It enables `zte` and `greedy`, disables `krylov`, and sets a proximity cutoff of 4.0. The initial density operator is `Lz` for `1H`. These settings describe the model; this script does not compare its output against measurement.

The wrapper sets `tmix=0.5`, offset 800, sweep `[1700 1700]`, `[512 512]` points, and `[2048 2048]` zero-fill sizes, and labels the plotted axes in ppm. It does not annotate units for the time, offset, or sweep literals, so those values are stated without assigning additional units.

The `liquid(spin_system,@roesy,parameters,'nmr')` call delegates the sequence to `@roesy` and supplies cosine and sine signal components. Each is apodised with a squared-cosine window in both dimensions; dimension-1 Fourier transforms are combined as `f1_cos-1i*f1_sin` and Fourier-transformed along dimension 2. The real spectrum is plotted with `plot_2d`. No pulse timings/phases, gradients, or receiver settings are exposed in this wrapper, so the sequence internals are not described.
