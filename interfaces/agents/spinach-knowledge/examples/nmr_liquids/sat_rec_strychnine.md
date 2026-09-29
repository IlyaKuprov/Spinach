# examples/nmr_liquids/sat_rec_strychnine.m

[Spinach source](https://github.com/IlyaKuprov/Spinach/blob/main/examples/nmr_liquids/sat_rec_strychnine.m)

This wrapper simulates a proton saturation-recovery experiment for strychnine. The source comment describes the example as being at 250 MHz; the implementation obtains a proton spin system from `strychnine({'1H'})` and sets `sys.magnet=5.9` (5.9 T by Spinach convention). It does not read a measured spectrum or an experimental acquisition file. Its comment estimates minutes of calculation time; this was not timed here.

The relaxation setup is Redfield with `inter.equilibrium='dibari'`, `inter.rlx_keep='kite'`, `tau_c={200e-12}` s (200 ps), and `temperature=298` (the wrapper gives no unit annotation). The basis is `sphten-liouv` / `IK-2`, with scalar-coupling connectivity and proximity level 1; the proximity cutoff is 5.0 and `greedy` parallelisation is enabled.

The single-channel acquisition uses `1H`, offset 1250, sweep 2500, and 4096 points, with ppm axis labels and an inverted axis. The wrapper sets a maximum delay of 0.5 and requests 10 delays; it does not state the delay list or explicitly label the units of these parameter literals. The call `liquid(spin_system,@sat_rec,parameters,'nmr')` delegates the saturation-recovery pulse sequence to `@sat_rec`. The wrapper does not specify pulse timings/phases, gradient use, or receiver settings, so none are inferred.

The returned FIDs are apodised exponentially with parameter 6, Fourier-transformed along dimension 1, shifted with `fftshift`, and plotted as their real part with `plot_1d`. This source produces a plotted simulated delay-series spectrum; it does not report fitted recovery constants or measured data.