# examples/esr_sol_pulsed/hpa_triplet.m

- Signature: `hpa_triplet()`
- Source: [`examples/esr_sol_pulsed/hpa_triplet.m`](https://github.com/IlyaKuprov/Spinach/blob/main/examples/esr_sol_pulsed/hpa_triplet.m)

## Physical aim and spin model

This hypothetical X-band pulse-acquire ESR calculation models the powder spectrum of a photogenerated pentacene triplet at 0.33 T. It contains one spin-1 triplet electron (`E3`), with an isotropic Zeeman tensor `diag([2,2,2])` and zero-field splitting `D = 1360.1 MHz`, `E = −47.2 MHz`; both values enter `zfs2mat` in Hz, with zero Euler angles. The full Zeeman Hilbert-space basis is used without approximation. No other spins or relaxation terms are specified.

## Pulse-acquire protocol and processing

For each orientation on `rep_2ang_6400pts_sph`, the initial triplet state is built by `zftrip` from the population weights `[0.56, 0.31, 0.13]`, the rotated ZFS tensor, and the field-dependent Zeeman tensor. The driver sets the detected state to `L+`, the hard-pulse operator to `Ly`, and the flip angle to `π/4`; `hp_acquire` applies the pulse and computes the free-induction signal. The carrier offset is zero, the sweep width is 4 GHz, and 128 points are acquired on a GHz lab-frame axis, inverted for plotting. The powder-averaged FID is apodised with `crisp`, zero-filled to 512 points, Fourier transformed, and the real spectrum is plotted.

The script creates a figure; it does not save a spectrum file. The model assumes isotropic `g=2` and is explicitly hypothetical. The source header estimates a calculation time of seconds.
