# examples/quantum_tech/transmon_frog.m

Source: [examples/quantum_tech/transmon_frog.m](https://github.com/IlyaKuprov/Spinach/blob/main/examples/quantum_tech/transmon_frog.m)

- Signature: `transmon_frog()`

## Model and target

FROG here means Frequency Robust Gate. The source sets up a three-level transmon model (`T3`), with rotating-frame frequency parameter `0.5e6` (0.5 MHz) and anharmonicity `-295.1e6` (−295.1 MHz). It uses the Zeeman-Liouville formalism without approximation and obtains the drift Hamiltonian in the cavity context. No electron spin, defect, EPR selection rule, or dissipative operator is specified.

The task encoded in this example is state transfer: the initial density operator is the normalised `BL1` state, and the normalised target is formed from the three-level state vector `[1; -1i; 0]/sqrt(2)`. The two controls are the transmon quadratures `Cx` and `Cy`. An initial two-quadrature pulse is built from five sine coefficients per channel: `Fa=[-0.6137 -0.0247 0.0742 0.0507 0.0149]` and `Fb=[-0.0106 0.0334 0.0579 0.0140 -0.0416]`.

## Optimisation and plots

The pulse has 224 slices over 112 ns, with slice duration `t_g/nsteps`. The optimiser samples five pulse-power levels, `2*pi*[15e6 16e6 17e6 18e6 19e6]`, and uses the SNSA penalty with weight 0.01, the Goodwin method, and a 50-iteration limit. The source asks for control, spectrogram, and robustness plots during the optimisation, then calls `fmaxnewton` with `grape_xy`. These settings describe an optimisation setup; the source does not report a final fidelity, convergence result, or measured gate performance.

The source cites the FROG preprint via [DOI 10.48550/arXiv.2511.22580](https://doi.org/10.48550/arXiv.2511.22580).
