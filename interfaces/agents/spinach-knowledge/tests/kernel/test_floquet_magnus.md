# tests/kernel/test_floquet_magnus.m

Exercises the production finite-harmonic recurrence in both AHT and zero-mean Van Vleck gauges. Exact cases include a static generator, a commuting waveform, and a degenerate eigenspace compared by its orthogonal projector rather than individual eigenvectors.

For a seeded noncommuting three-dimensional matrix family, independent ode113 integration is tightened before comparison. The test checks non-stroboscopic error powers at two drive frequencies for effective powers zero through three, Hermiticity, unitarity, kick gauge conditions, and covariance under a joint Fourier-phase/time-origin shift. A complex-signal and full-propagator comparison detects omission of endpoint kicks. Negative modulation frequency and malformed inputs are covered. The random stream is restored on exit.

These are finite test cases, not a proof of general convergence or long-time error bounds.
