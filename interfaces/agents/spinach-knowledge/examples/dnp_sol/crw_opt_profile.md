# cRW-OPT static powder transfer

`crw_opt_profile()` uses Goodwin et al.'s waveform ([DOI](https://doi.org/10.1016/j.jmro.2026.100218)): 30 pieces of 5 ns (150 ns per block). The initial electron `Lx` is transferred to the proton, detected with `-Lz` without peak normalisation. The model is one electron and one proton at 14.8 MHz proton Larmor frequency with tilted axial hyperfine coupling and no relaxation.

The on-resonance build-up samples the 33 complete blocks in the original 5 μs window. The bandwidth plot uses 13 blocks at each signed electron detuning from −60 to +60 MHz in 2 MHz increments, with 100 powder nodes for the bandwidth and 400 for the more orientation-sensitive build-up. The reduced grids are for demonstration, not a claim of full-grid reproduction.
