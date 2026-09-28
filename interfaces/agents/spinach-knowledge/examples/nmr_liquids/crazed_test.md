# examples/nmr_liquids/crazed_test.m

- Signature: `crazed_test()`

## Purpose

Long-range intermolecular coherences predicted by Warren and co-workers. Calculation time: seconds. Reference: [doi:10.1126/science.8266096](http://dx.doi.org/10.1126/science.8266096).

## Physical / mathematical content

This example models a four-proton system with specified positions and temperature to demonstrate the long-range intermolecular coherence phenomenon associated with the CRAZED experiment. It prepares a thermal-equilibrium state and simulates the two-dimensional CRAZED signal.

## Numerical / algorithmic content

The code uses the complete `sphten-liouv` basis, sets the field to 6.0 and specifies a 90-degree angle, offset 1300, sweep 5000, 512 points and 2048 zero-fill points on each axis. The initial state is constructed in the lab frame from the Hamiltonian and equilibrium routine at orientation [0 0 0]. A cosine window is applied before a shifted 2D FFT; plotting uses the spectrum magnitude and the positive display mode.

## Implementation structure

The function defines four proton sites, shifts, coordinates and temperature 100, constructs the Spinach system and complete basis, prepares `rho0`, and calls `crystal(...,@crazed,...,'nmr')`. It then apodises and transforms the FID for plotting.
