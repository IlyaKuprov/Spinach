# fme_bloch_siegert.m

- Source: `examples/fundamentals/fme_bloch_siegert.m`
- Entry point: `fme_bloch_siegert()`; no inputs or outputs
- Reference: Bock, Ernst, and Uhrig, DOI:10.5194/mr-2026-10, Section 4.1, Eq. (36), Fig. 4

A single proton in Hilbert space demonstrates counterrotating-field micromotion in the resonant rotating frame. The carrier is 10 kHz, the nutation frequency 2 kHz, and the Fourier fundamental 20 kHz. Spinach operators construct the three signed harmonics. The leading longitudinal Bloch-Siegert coefficient is printed in Hz.

The order-2 Van Vleck generator and kicks are compared with direct propagator integration by ode113, with relative/absolute tolerances 1e-10/1e-12. The initial state and observable are Lz; the signal is normalised by trace(Lz^2). Complete evolution uses both endpoint kicks. A separate curve deliberately omits both kicks. Plots display the signal and errors over 40 field periods, marking period boundaries. This is an illustrative time-domain calculation, not a reproduction of the paper's Lx-detected spectra. Printed maximum errors are observations, not assertions of universal superiority or convergence.
