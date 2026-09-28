# examples/extremes/phosphorus_cluster.m

- Signature: `phosphorus_cluster()`

## Purpose

Phosphorus system simulation for Gerhard Hagele using brute-force Liouville-space time propagation. WARNING: needs 32+ CPU cores, 128+ GB of RAM, and a strong FP64-capable Nvidia GPU. Run time on the above: hours.

## Physical / mathematical content

- The 34-spin system contains seven `31P` spins and 27 `1H` spins at magnet induction `7.04`. Phosphorus chemical shifts are `-99.6`, `-0.33` (spins 2–4), and `-156.7` (spins 5–7); proton shifts are `0.21` (spins 8–34).
- The specified phosphorus–phosphorus J-couplings are `-323.22`, `46.18`, `-354.82`, `25.79`, `-9.04`, `-16.63`, and `-214.10` across the spin pairs assigned in the source. Each of spins 2, 3, and 4 couples at `4.0` to its respective nine-proton group: 8–16, 17–25, or 26–34.

## Numerical / algorithmic content

- The basis uses `sphten-liouv` formalism, `IK-2` approximation, `scalar_couplings` connectivity, proximity level `1`, longitudinal `1H` spins, and projection `1`. Three `S3` symmetry groups act on spins `[8 9 10]`, `[17 18 19]`, and `[26 27 28]`. Greedy parallelisation is enabled.
- The `31P` acquisition uses `L+` for both the initial state and detection coil, with no decoupling, offset `-10000`, sweep `30000`, `16536` points, zero fill to `65536`, an axis in `Hz`, and axis inversion `1`.

## Implementation structure

- The function creates the spin system and basis, runs `liquid` with `@acquire` in `nmr` mode, applies exponential apodisation with parameter `5`, computes `real(fftshift(fft(fid,parameters.zerofill)))`, and plots the resulting one-dimensional spectrum.