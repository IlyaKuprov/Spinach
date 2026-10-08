# examples/nmr_liquids/roesy_strychnine.m

[Spinach source](https://github.com/IlyaKuprov/Spinach/blob/main/examples/nmr_liquids/roesy_strychnine.m)

This example simulates a proton ROESY spectrum for a strychnine spin system. It obtains the system through `strychnine({'1H'})`, rather than loading a measured spectrum or molecular-structure file in this wrapper. It sets `sys.magnet=5.9` (5.9 T by Spinach convention) and uses a `sphten-liouv` basis with the `IK-2` approximation, scalar-coupling connectivity, and proximity level 3. The relaxation model is Redfield, with zero equilibrium, secular relaxation terms, and `tau_c={200e-12}` (200 ps). The wrapper enables `zte` and `greedy`, disables `krylov`, and sets the proximity cutoff to 4.0.

The initial state is `Lz` for `1H`. Sequence parameters are `tmix=0.5`, offset 1200, sweeps `[2500 2500]`, `[512 512]` points, and `[2048 2048]` zero-fill sizes. It requests ppm axes. The wrapper supplies no explicit units for the `tmix`, offset, and sweep literals; they are reported as passed, not reinterpreted here.

The call `liquid(spin_system,@roesy,parameters,'nmr')` delegates the ROESY sequence to `@roesy` and returns cosine and sine signal components. Both receive squared-cosine apodisation in both dimensions. The code Fourier-transforms the components along dimension 1, selects the imaginary part of the cosine signal and real part of the sine signal, combines these as a States quadrature signal, then Fourier-transforms dimension 2. It plots the real spectrum with `plot_2d`. This wrapper gives no pulse list, gradient scheme, or receiver-phase settings; the pulse-program internals are not specified here. The header estimates runtime as minutes; that is a comment, not a runtime measured for this draft.
