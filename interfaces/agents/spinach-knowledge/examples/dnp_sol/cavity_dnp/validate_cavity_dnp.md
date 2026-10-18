# validate_cavity_dnp

`validation=validate_cavity_dnp()` executes the native `sphten-liouv` example at cavity cutoffs 3, 4, and 5, repeats propagation on a half-step grid, and checks trace, pulse swaps, population bounds, and unoccupied Fock edge. It also tests a zero-exchange control, lossless excitation conservation and resonant Rabi oscillations, and the exact resolved damped Jaynes–Cummings photon transient.

All simulated controls retain the same formalism and use native Hamiltonian, state/operator, and evolution APIs. Analytical expressions are independent reference curves only. Success prints `CAVITY_DNP_VALIDATION_SUCCESS`; the function writes no files. These tests validate the documented four-level analogue, not the omitted low-field Bi donor or ensemble physics.
