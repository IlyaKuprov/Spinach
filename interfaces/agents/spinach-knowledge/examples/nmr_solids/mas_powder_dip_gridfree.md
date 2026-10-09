# examples/nmr_solids/mas_powder_dip_gridfree.m

Source: [examples/nmr_solids/mas_powder_dip_gridfree.m](https://github.com/IlyaKuprov/Spinach/blob/main/examples/nmr_solids/mas_powder_dip_gridfree.m)

## Model

This source describes a powder MAS spectrum for a pair of dipole-coupled protons using a grid-free Fokker-Planck equation. It sets a 14.1 T field, two `1H` isotropic Zeeman scalar values (5.0 and -2.0), and coordinates [0, 0, 0] and [0, 3.9, 0.1]. The coordinates define the dipolar geometry on system construction; the code supplies no separate dipolar-coupling constant. Coordinate and Zeeman scalar units are not given in the file. There is no quadrupolar spin or NQI, and no explicit RF pulse or decoupling field.

The basis uses spherical-tensor Liouville space, no approximation, and projection +1. The simulation settings include MAS axis [1, 1, 1], rotor rate 1000 Hz, maximum rank 15, and 512 acquisition points. The algorithm is specifically selected by `gridfree(spin_system,@acquire,parameters,'nmr')`; this is distinct from the `floquet` and `singlerot` calls in the paired examples.

## Signal and display

The initial state and receiver coil are both the proton `L+` state. The sweep setting is 2e4, offset is 0, zero filling is 4096, and the plotted axis is configured in ppm with inversion enabled. After `gridfree` produces the FID, the code applies exponential apodisation parameter 6, Fourier transforms using `fftshift`, and plots the real spectrum through `plot_1d`. The plot is a simulation output, not an experimental measurement reported by the source.
