# examples/esr_sol_pulsed/holeburn_gd_dota_powder.m

- MATLAB implementation: [examples/esr_sol_pulsed/holeburn_gd_dota_powder.m](https://github.com/IlyaKuprov/Spinach/blob/main/examples/esr_sol_pulsed/holeburn_gd_dota_powder.m)

[MATLAB example](../../../../../examples/esr_sol_pulsed/holeburn_gd_dota_powder.m) · [holeburn sequence helper](../../../../../experiments/holeburn.m)

- Signature: `holeburn_gd_dota_powder()`

## Aim and model

This example models soft-pulse spectral hole burning for a powder of Gd(III) centres. It samples a zero-field-splitting (ZFS) distribution using `zfs_sampling(30,5,1e-2)`; the source says its statistical parameters come from Figure 5 of Raitsimring et al. The run samples each returned `D,E` pair with weight `W`, and applies a numerical powder grid and a numerical second-order rotating-frame transformation.

For each ZFS sample, the system is `E8` (the electron-spin label used for Gd(III)) at 3.5 T, with isotropic Zeeman scalar 2.002319 and a ZFS matrix formed as `0.56e9*zfs2mat(D(n),E(n),0,0,0)`. The basis is spherical-tensor Liouville space, exact (`approximation='none'`), with projections `-3:3`; trajectory-level SSR algorithms are disabled. The factor and arguments are recorded as coded; the source does not state separate units for the ZFS sampler outputs.

## Pulse and acquisition protocol

The sequence helper `holeburn` applies the shaped soft pulse using the Fokker–Planck formalism, then an ideal hard `pi/2` observation pulse and signal acquisition. The soft-pulse rank is 2, phase `-pi/2` rad, carrier `-0.5e9` Hz, duration `50e-9` s, and propagator method `expm`. Two otherwise matched simulations are compared: A sets the soft-pulse power to zero; B sets it to `2*pi*1e7` rad/s. Both use the same pulse frequency, phase and duration.

The initial state is `Lz` on `E8`, the detection operator is `L+`, and no spins are decoupled. Acquisition uses zero offset (the helper does not annotate offset units), sweep width `0.8e10` Hz, 512 points, a 2048-point zero-fill, the `rep_2ang_400pts_sph` orientation grid, second-order rotating frames for `E8`, no derivative, and a non-inverted axis. The axis display unit is GHz. Sweep width is in Hz as defined by the `acquire` helper; pulse parameter units follow `holeburn`'s parameter documentation.

## Observable and output

For each ZFS sample, the two FIDs receive exponential apodisation with parameter 10, are Fourier transformed, and are accumulated with that sample's weight. The real accumulated spectrum for A is plotted in red and B in blue, with the plot refreshed during the weighted ZFS loop. The example plots the spectra and does not write a data or figure file.

The source notes that non-central Gd(III) transition holes are very shallow and estimates a calculation time of minutes. The cited distribution source is Raitsimring et al., *Applied Magnetic Resonance* 28, 281–295 (2005), Figure 5 ([DOI](https://doi.org/10.1007/BF03166762)).
