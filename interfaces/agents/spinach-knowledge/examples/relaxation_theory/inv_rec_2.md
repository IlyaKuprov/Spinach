# examples/relaxation_theory/inv_rec_2.m

- MATLAB implementation: [examples/relaxation_theory/inv_rec_2.m](https://github.com/IlyaKuprov/Spinach/blob/main/examples/relaxation_theory/inv_rec_2.m)

Source: [examples/relaxation_theory/inv_rec_2.m](https://github.com/IlyaKuprov/Spinach/blob/main/examples/relaxation_theory/inv_rec_2.m)

## Purpose and molecular model

Simulates proton inversion-recovery spectra for the strychnine spin system at six recovery delays; the source estimates calculation time in minutes. It obtains the proton system from strychnine({'1H'}) and sets the field to 14.1 T. It is a simulation workflow, not an experimental result or a claimed fit.

## Relaxation and basis

The relaxation mechanism is Redfield with tau_c={200e-12} s (200 ps), dibari equilibrium, kite retention (rlx_keep='kite'), and temperature value 298 (the source does not state a temperature unit). The basis uses sphten-liouv, IK-2, scalar-coupling connectivity, and proximity level 3. The script sets the proximity cutoff to 4.0 without stating its unit and disables Krylov propagation.

## Pulse sequence and acquisition

The recovery delays are 0.01, 0.1, 0.5, 1, 5, and 10 s. For each delay, the script starts from isotropic thermal equilibrium, applies a pi rotation about Ly, evolves under the rotating-frame Hamiltonian plus Redfield relaxation for that delay, then applies a pi/2 Ly read pulse. It acquires the proton signal with an L+ detection state. Acquisition settings are sweep 6500, 8192 points, zero-fill to 65536, proton channel, ppm axis, and an explicitly labelled offset of 2800 Hz. Each FID receives exponential apodisation with parameter 5, is Fourier transformed, and is plotted in its own panel of a 2-by-3 figure labelled by recovery delay. The plot therefore compares six simulated spectra across the specified delays; the source does not report a measured spectrum or rate benchmark.
