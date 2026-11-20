# Powder stimulated-echo ESEEM (`eseem_3p_powder`)

Run `eseem_3p_powder()` to compare the in-phase stimulated echoes of a model electron–proton pair for ideal and 25 ns rectangular 90° pulses. The transverse electron preparation and receiver are both `Lx`, and the outer delays are 104 ns. The storage delay is sampled every 4 ns.

The `esr` assumption set retains the proton Larmor precession and the electron-secular part of the tilted axial hyperfine tensor. `powder` integrates the complex echo over static orientations at each carrier offset; the example then takes the arithmetic mean over a uniform ±15 MHz electron-offset distribution. Keep this averaging ahead of any magnitude or Fourier processing. The example uses a native spherical quadrature rather than reproducing a particular publication grid.
