# examples/extremes/phosphorus_cluster.m

- MATLAB implementation: [examples/extremes/phosphorus_cluster.m](https://github.com/IlyaKuprov/Spinach/blob/main/examples/extremes/phosphorus_cluster.m)

- Entry point: `phosphorus_cluster()`.
- Source: [`examples/extremes/phosphorus_cluster.m`](../../../../../examples/extremes/phosphorus_cluster.m).

## System and scientific intent

This is a 31P NMR simulation for a phosphorus system attributed in the source comments to Gerhard Hagele. The isotope list contains 7 `31P` and 27 `1H` (34 spins total). The script sets the magnetic-induction parameter to `7.04` and explicitly assigns the Zeeman and scalar-coupling tables; their entries are not given units in the source comments. For example, the first phosphorus spin is coupled to spins 2–4 with value `-323.22` and to spins 5–7 with `46.18`; the remaining P–P and P–H couplings are also assigned directly in the source.

The basis uses `sphten-liouv`, `IK-2` approximation, scalar-coupling connectivity and proximity level `1`, with longitudinal `1H` terms and projection `1`. Three `S3` groups are specified for proton triplets `[8 9 10]`, `[17 18 19]` and `[26 27 28]`. Greedy parallelisation is enabled; the commented GPU option is not active in the `sys.enable` assignment.

## Acquisition and plotted observable

The sequence observes `31P` with `L+` initial state and receiver, no decoupling, offset `-10000` and sweep `30000`. It uses `16536` points, zero-fill to `65536`, output axis units `Hz`, and inverted axis. The source does not separately annotate offset/sweep units. No explicit RF-pulse table is defined: the code calls `liquid(spin_system,@acquire,parameters,'nmr')`.

Exponential apodisation uses parameter `5`; the FID is Fourier-transformed with zero-fill length, made real, and plotted as a one-dimensional spectrum. The source contains no numerical peak positions, intensities, or comparison against experiment.

## Practical limit

The source describes brute-force Liouville-space time propagation and warns of 32+ CPU cores, 128+ GB RAM, a strong FP64-capable Nvidia GPU, and a run time of hours.

Zero track elimination is explicitly enabled with `zte` in `sys.enable`.
