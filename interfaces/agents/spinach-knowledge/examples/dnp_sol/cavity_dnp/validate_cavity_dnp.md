# validate_cavity_dnp

`validation=validate_cavity_dnp(reference)` validates the actual E/209Bi/Cn example, with an optional previously computed default reference. It compares Fock cutoffs 3, 4, and 5, increases the pulse Fokker–Planck rank, doubles the RF saturation length, and checks half-size and non-dyadic laboratory-frame reset timesteps.

Independent native controls omit the RF pulse or the spin–cavity exchange, test lossless conservation of total spin projection plus photon number, resolve the donor-generated photon transient, and compare native one-photon loss with the exact exponential. Trace, observable bounds, saturation residual, and Fock-edge occupation are checked. Every simulated trajectory stays in `sphten-liouv`; no density-matrix conversions, custom generators, or effective three-level dynamics are used. Success prints `CAVITY_DNP_VALIDATION_SUCCESS`.
