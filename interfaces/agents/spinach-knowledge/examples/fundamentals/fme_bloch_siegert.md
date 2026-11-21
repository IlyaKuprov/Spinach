# fme_bloch_siegert.m

- Source: `examples/fundamentals/fme_bloch_siegert.m`
- Entry point: `fme_bloch_siegert()`; no inputs or outputs
- Reference: Bock, Ernst, and Uhrig, DOI:10.5194/mr-2026-10, Section 4.1, Eq. (36), Fig. 4

A single proton in Hilbert space demonstrates counterrotating-field micromotion in the resonant rotating frame. The carrier is 10 kHz, the nutation frequency 2 kHz, and the Fourier fundamental 20 kHz. Spinach operators construct the three signed harmonics. The leading longitudinal Bloch-Siegert coefficient is printed in Hz.

The order-2 Van Vleck generator and kicks are compared with native fourth-order three-point `step` quadrature of the physical Hamiltonian. Eight and sixteen substeps per plotted interval are compared; the example requires their normalised signals to agree within 1e-8. `state` and `coil_state` construct the initial Lz state and detector, and `hdot` evaluates the observable and its initial normalisation.

Effective-generator trajectories use `evolution`, and `step` applies each dimensionless endpoint kick. The deliberately kick-free curve uses native observable evolution. Plots retain the same 40 field periods and period-boundary markers. This is an illustrative time-domain calculation, not a reproduction of the paper's Lx-detected spectra. Printed approximation errors are observations, not assertions of universal superiority or convergence.
