# validate_cavity_dnp

`validation=validate_cavity_dnp(reference)` checks a result from `cavity_dnp_minimal(3,2,300)`. A no-argument call first runs that reference calculation.

Checks cover cavity dimensions 3, 4, and 5; pulse phase-grid rank; saturation-sweep count; and reset time-step refinements by factors of 2 and 5/4. Physical controls omit the selective RF pulse or spin-cavity exchange, test lossless conservation of donor projection plus photon number, resolve the photon transient, and compare one-photon loss with the exact exponential. Trace, observable bounds, saturation residual, and Fock-edge occupation are also checked.

All trajectories remain in `sphten-liouv`. The function returns the measured errors and diagnostics and prints `CAVITY_DNP_VALIDATION_SUCCESS` after all assertions pass. It writes no files or figures.
