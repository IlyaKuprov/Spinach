# examples/quantum_tech/transmon_frog.m

- Signature: `transmon_frog()`

## Purpose

A Frequency Robust Gate (FROG) for a single transmon, using ensemble GRAPE over a distribution of control powers with an excess-amplitude penalty. Calculation time: minutes. The source cites https://doi.org/10.48550/arXiv.2511.22580.

## Model and parameters

The T3 transmon is represented in a rotating frame with frequency 0.5 MHz and anharmonicity -295.1 MHz, using the Zeeman-Liouville formalism without approximation. The gate duration is 112 ns with 224 control steps. Its target state is the normalized three-level superposition `[1; -1i; 0]/sqrt(2)`.

## Optimization

The initial two-quadrature pulse is built from five sine coefficients per channel. GRAPE optimizes it over power levels 15–19 MHz, with SNSA penalty weight 0.01, the Goodwin method, and a 50-iteration limit. The source plots controls, a spectrogram, and robustness during optimization, and runs `fmaxnewton` with `grape_xy`.
