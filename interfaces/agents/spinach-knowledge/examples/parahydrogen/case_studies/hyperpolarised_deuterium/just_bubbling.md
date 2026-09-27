# examples/parahydrogen/case_studies/hyperpolarised_deuterium/just_bubbling.m

- Signature: `just_bubbling()`

## Purpose

Simulates the evolution of state populations during ortho-deuterium bubbling in the presence of a parahydrogenation catalyst, followed by free evolution. No pulses are applied. A paper link is to follow.

## Spin system and parameters

- Four `2H` spins with experimental chemical shifts `{4.55, 4.55, -13.5, -16.5}` and J-couplings `J(1,2) = 12.0`, `J(3,4) = 0.24`.
- NQI tensors for spins 3 and 4, from a DFT calculation, are respectively `1e3*[108.2 0.1 28.1; 0.1 -55.6 4.5; 28.1 4.5 -52.6]` and `1e3*[-55.2 8.1 -14.5; 8.1 -6.3 -73.9; -14.5 -73.9 61.5]`.
- DFT Cartesian coordinates are unspecified for the D₂ pair; spins 3 and 4 have coordinates `[-1.98 0.45 -0.55]` and `[-0.25 1.33 -1.65]`.
- Chemical-exchange parts are `{[1 2], [3 4]}`, with rate matrix `[-1 5000; 1 -5000]` and initial concentrations `[1 0]`. The magnetic field is `7.05`.
- The simulation uses the `sphten-liouv` formalism with approximation `none`. Relaxation settings are `{'redfield','t1_t2'}`, equilibrium `zero`, secular terms retained, correlation times `{1e-12, 400e-12}`, R₁ rates `{0.04 0.04 0 0}`, and R₂ rates `{8.00 8.00 0 0}`.

## Evolution and output

The initial state is `unit_state(spin_system)`. The singlet, triplet, and quintet states of spins 1 and 2 are obtained with `deut_pair(spin_system,1,2)`; their traceless components are used for detection. After applying the `nmr` assumptions, the free-evolution Liouvillian is `LF=H+1i*R+1i*KF`. Bubbling uses `LB=H+1i*R+1i*KB`, where `KB` is produced by `magpump` from a traceless `S+Q{1}+Q{2}+Q{3}+Q{4}+Q{5}` state at the provisional rate `1e-1`.

Bubbling runs for 7 seconds with a `0.007` time step and `1000` steps; free evolution then runs for 30 seconds with a `0.03` time step and `1000` steps. The script plots the real coefficients of `S`, `T_{\pm1}`, `T_{0}`, `Q_{\pm2}`, `Q_{\pm1}`, and `Q_{0}` against time.