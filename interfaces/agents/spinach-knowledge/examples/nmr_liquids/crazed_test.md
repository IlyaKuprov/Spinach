# examples/nmr_liquids/crazed_test.m

- MATLAB implementation: [examples/nmr_liquids/crazed_test.m](https://github.com/IlyaKuprov/Spinach/blob/main/examples/nmr_liquids/crazed_test.m)

- Signature: `crazed_test()`

## Purpose

A four-spin CRAZED simulation illustrating long-range intermolecular coherences discussed by Warren and co-workers. The source cites [doi:10.1126/science.8266096](http://dx.doi.org/10.1126/science.8266096) and estimates a calculation time of seconds.

## Spin system and preparation

The model contains four `1H` sites at a 6.0 T field, with chemical shifts 2.0, 2.0, 8.0, and 8.0 ppm. Four position vectors are entered as [0,0,0], [1/2,1/2,-1/sqrt(2)] times 1e2, [1,0,0] times 1e2, and [1/2,-1/2,-1/sqrt(2)] times 1e2. The example sets the spin-system temperature to 100 K; it does not state coordinate units alongside these vectors. It uses the complete Liouville-space basis (no basis approximation).

The initial state is constructed with `equilibrium` from the Hamiltonian and Q operators in the lab frame, using [0,0,0] as the final argument. The sequence parameters set a pi/2 pulse angle, 1300 Hz offset, 5000 Hz sweep width, and orientation [0,0,0].

## CRAZED acquisition and processing

The crystal simulation samples a 512 by 512 FID, then zero-fills to 2048 by 2048 for the 2D FFT. Cosine apodisation is applied along both dimensions. The plotted data are the spectrum magnitude, with positive contours; the axes are not assigned explicit units in this function.

## Interpretation and scope

The code prepares one specified four-spin model and orientation and calls the CRAZED sequence through the crystal-simulation path. Its displayed magnitude spectrum is a simulation output, not a measured spectrum or a powder-averaged result. The example does not supply a separate experimental comparison, and coordinate units are not identified in the source.
