# examples/esr_sol_pulsed/hpa_triplet.m

- Signature: `hpa_triplet()`

## Purpose

Simulates a hypothetical powder-averaged X-band pulse-acquire ESR spectrum of a photogenerated pentacene triplet at 0.33 T. Calculation time: seconds.

## Physical model

The model uses one spin-1 triplet (`E3`) with an isotropic g tensor (g = 2) and zero-field-splitting parameters D = 1360.1 MHz and E = −47.2 MHz. The ZFS principal axes are supplied with zero Euler angles. The initial triplet populations are [0.56, 0.31, 0.13]; the `zftrip` state is rotated with the powder orientation and includes the field-dependent Zeeman term.

## Simulation and processing

The calculation uses the full Zeeman Hilbert-space basis and disables trajectory-level SSR. It applies an `Ly` pulse of angle π/4, then powder-averages the pulse-acquire ESR response on the `rep_2ang_6400pts_sph` grid. Acquisition uses a 4 GHz sweep, 128 points, zero offset, and a GHz lab-frame axis; the axis is inverted. The FID is apodised with `crisp`, zero-filled to 512 points, Fourier transformed, and the real spectrum is plotted.
