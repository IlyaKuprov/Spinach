# cce_noise_spectrum

Syntax: [answer,diagnostics]=cce_noise_spectrum().

A mixed carbon-13 nuclear spin in a 0.01 T field supplies the ordered,
connected transverse fluctuation correlation. The example defines ordinary
Spinach system tensors and a named Lx fluctuation descriptor, then calls
clexpn(spin_system,@noise_correlation,parameters,'labframe') once. The
spectrum is returned in diagnostics.spectrum; plots retain the same time
and angular-frequency grids, Hann window, trapezoidal quadrature, and
stationary two-sided Fourier convention. No hand-built operator or state is
required. The independent cosine reference and its numerical error are
checked in test_cce_examples, not in the example.
