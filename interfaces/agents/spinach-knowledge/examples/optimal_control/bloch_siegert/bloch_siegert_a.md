# examples/optimal_control/bloch_siegert/bloch_siegert_a.m

- Signature: `bloch_siegert_a()`
- Source: [examples/optimal_control/bloch_siegert/bloch_siegert_a.m](https://github.com/IlyaKuprov/Spinach/blob/main/examples/optimal_control/bloch_siegert/bloch_siegert_a.m)

## Objective and spin model

Estimated calculation time: minutes. The example optimises a 90-degree state transfer, `Lz -> Lx`, for a single on-resonance `1H` spin at `sys.magnet=14.1`. The scalar Zeeman interaction is set to zero, and the calculation uses the exact `sphten-liouv` basis. The initial and target states are constructed and normalised separately. The drift is the NMR-assumption Hamiltonian; the two control operators are `Lx` and `Ly` on the proton channel.

## Waveform parameterisation and comparison

For each of 20 control levels, the example scales the absolute proton Larmor angular frequency by evenly spaced relative powers from 0.001 to 1.0. Power levels are in rad/s. Each pulse has 50 equal-duration slices, with slice duration `(pi/100)/pwr_level`; the common initial guess is a 2-by-50 Gaussian array scaled by 1/10. The optimiser is L-BFGS through `fmaxnewton`, with at most 500 iterations and `tol_x=1e-4`.

At each power it constructs settings with Bloch-Siegert (BSS) correction disabled and enabled, then optimises one `grape_xy` waveform for each setting. Both waveforms are evaluated with the BSS-enabled settings via `ensemble`. The plotted quantity is terminal infidelity, `1-fidelity`, against relative control power, with a logarithmic vertical axis. Both designs are therefore evaluated in the presence of BSS physics.
