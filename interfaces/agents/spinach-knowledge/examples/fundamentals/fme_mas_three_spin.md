# fme_mas_three_spin.m

- Source: `examples/fundamentals/fme_mas_three_spin.m`
- Entry point: `fme_mas_three_spin()`; no inputs or outputs
- Reference: Bock, Ernst, and Uhrig, DOI:10.5194/mr-2026-10, Sections 4.2-4.3, Eq. (39), Figs. 9-11

3 protons at 9.4 T with isotropic shifts 0, 2, and 5 ppm and point-dipolar couplings undergo 25 kHz MAS. Coordinates in Angstrom are [0 0 0], [0 0 2.5], and [2.2 0.6 1.0]. One crystallite has laboratory Euler angles [0.3 0.7 0.2] radians. These are illustrative parameters, not the paper's spin system or powder average. Positive native MAS rate is represented by angular frequency `-2*pi*rate` in the Fourier reconstruction.

Native rotor_stack returns the five physical Hamiltonians required to reconstruct the rank-two signed Fourier coefficients using MATLAB fft and fftshift. This reuses Spinach orientation and interaction handling. The existing drifts exporter with singlerot returns a spatially enlarged Fokker-Planck generator; rotor_stack supplies the physical blocks without constructing and subsequently unpacking that enlarged matrix.

At inverse-frequency order 2, complete AHT and Van Vleck propagators, each with both endpoint kicks, are compared against native three-point step quadrature of the Fourier Hamiltonian. Eight and sixteen substeps per plotted interval must give normalised signals within 1e-8. Native state and coil_state construct total Lx preparation and detection; hdot supplies the observable and initial normalisation. Effective trajectories use evolution, and step applies the dimensionless endpoint kicks. Plots show 20 rotor periods with 40 intervals per rotor and explicitly mark rotor-synchronised samples. Both gauges retain their periodic kicks at within-rotor points as well as at stroboscopic times. The example reports errors rather than asserting convergence or reproducing spectral figures.
