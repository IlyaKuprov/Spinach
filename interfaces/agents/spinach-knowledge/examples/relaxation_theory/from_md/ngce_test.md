# examples/relaxation_theory/from_md/ngce_test.m

- MATLAB implementation: [examples/relaxation_theory/from_md/ngce_test.m](https://github.com/IlyaKuprov/Spinach/blob/main/examples/relaxation_theory/from_md/ngce_test.m)

- Signature: `ngce_test()`
- Run from MATLAB with no arguments. It prints the analytical/numerical relaxation matrices and opens a comparison figure; it returns no MATLAB value.
- Calculation time stated by the source: minutes.

## Purpose

This is a numerical check of Spinach's `ngce` numerical-integral route to a Redfield relaxation superoperator against the analytical Redfield matrix for isotropic rotational diffusion. It is a stochastic simulation-based comparison, not a deterministic unit test with a pass/fail threshold.

## Model and numerical settings

The source uses a two-spin `{'1H','13C'}` system at `sys.magnet=9.4`, with coordinates `[0 0 0]` and `[0 0 1.02]`. It requests Redfield relaxation, zero equilibrium, lab-frame retention, `inter.tau_c={1.0e-10}`, and the full `sphten-liouv` basis (`bas.approximation='none'`). The file does not annotate physical units for these literal field, coordinate, and time values; use the conventions expected by Spinach rather than inferring units from the example.

## Numerical comparison

After constructing the spin system and basis, the example takes the real part of the analytical matrix `relaxation(spin_system)`. It obtains lab-frame `H0` and interaction components `Q`, then samples `rwalk(100000,tau_c,tau_c/25)` (so `dt=tau_c/25`). Each Euler-angle sample is converted with `orientation(Q,...)` into a Hamiltonian-trajectory cell entry; this loop is run with `parfor`. The trajectory and matrices are passed to `ngce(spin_system,H0,H1,dt,tau_c,1e-3)`; the meaning/units of the final `1e-3` argument are not explained in this example.

The command window receives the full analytical matrix, numerical matrix, and `dR_gce` matrix labelled as the standard deviation of the mean. The plotted points compare the diagonal rates (numerical `gce_rates` on x, analytical `red_rates` on y) against a unity-slope reference; horizontal error bars are `2*diag(dR_gce)`. The source labels these as 95% confidence bands. Because the walk is random and this function sets no random seed, repeated runs need not give identical numerical estimates. No assertion or acceptance threshold is applied.
