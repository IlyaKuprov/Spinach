# examples/optimal_control/features_dt_opt.m

[Stable source link](https://github.com/IlyaKuprov/Spinach/blob/main/examples/optimal_control/features_dt_opt.m) · [Pulse-sequence reference](https://doi.org/10.1016/0022-2364(83)90133-6)

## Calculation

This example changes only the six slice durations of a composite 13C inversion pulse: the Cartesian control amplitudes and phases remain fixed. The starting sequence, also given in Fig. 3 of the cited paper, is `270(-x)360(x)90(y)270(-y)360(y)90(x)`. The simulated ensemble contains 100 non-interacting 13C spins with chemical shifts evenly spanning −166 to +166 ppm at 14.1 T (about ±25 kHz as described in the source). The `sphten-liouv` `IK-2` basis uses proximity level 1 and scalar-coupling connectivity; there are no couplings between these spins. The normalised initial state is 13C `Lz` and is also supplied as the target; the source comments identify minimising the objective as inversion.

The six-slice Cartesian waveform is `2π×25 kHz` times x-coefficients [−1, 1, 0, 0, 0, 1] and y-coefficients [0, 0, 1, −1, 1, 0]. Initial durations are [30, 40, 10, 30, 40, 10] μs, totaling 160 μs. The `tgrape` objective supplies the gradient to MATLAB `fmincon`; an L-BFGS Hessian approximation is requested, all durations have lower bound zero, and an equality constraint fixes their sum at 160 μs. Thus the optimiser redistributes time among the six slices without changing the waveform amplitude or overall duration.

## Simulated observables

The source propagates both the initial and optimised durations with `step`. For each, it applies a 90° `Ly` readout rotation, simulates 13C acquisition with `liquid`/`acquire`, applies Gaussian apodisation (parameter 10), zero-fills to 16384 points, and forms the shifted FFT on an inverted Hz axis. The script plots initial and optimised spectra side by side and prints the duration vectors in microseconds. The DOI and the source comment's “slightly better” design motivation are retained, but the file contains no numerical improvement or experimental validation; none is asserted here.
