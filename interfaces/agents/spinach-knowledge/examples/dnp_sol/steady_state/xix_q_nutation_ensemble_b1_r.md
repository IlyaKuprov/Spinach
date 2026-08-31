# examples/dnp_sol/steady_state/xix_q_nutation_ensemble_b1_r.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/examples/dnp_sol/steady_state/xix_q_nutation_ensemble_b1_r.m`
- Signature: `xix_q_nutation_ensemble_b1_r()`
- Total lines: 134

## Purpose

Simulation of nutation frequency dependence of XiX DNP field profiles in the steady state with electron-proton distance and electron Rabi frequency ensembles. Calculation time: minutes

## Physical / mathematical content

- Solid-state DNP examples. These files model microwave-driven electron-nuclear polarisation transfer mechanisms such as the solid effect, cross effect, NOVEL, XiX, TOP, BEAM, and TPPM variants. The mathematics combines driven spin dynamics, relaxation, powder/MAS averaging, and steady-state or transient propagation.

## Numerical / algorithmic content

- The file is built around the standard Spinach workflow: create the spin system, choose a basis or context, assemble operators/superoperators, then propagate or analyse the resulting dynamics.
- The file also defines local helper function(s): `xix_field_profile_b1_r()`. This usually means the public entry point is supported by tightly coupled validation or helper logic kept private to the file.

## Code-derived implementation details

### Comment-guided execution stages

- Lines 13-14: Nutation frequencies, Hz; implemented by `nu=1e6*[6.8 9.6 13.5 17.5 25 36]`.
- Lines 16-17: Shot repetition times, seconds; implemented by `srt=1e-3*[0.051 0.051 0.102 0.153 0.153 0.306]`.
- Lines 19-20: Get the figure started; implemented by `kfigure(); hold on; kgrid`.
- Lines 28-29: Plot the curves; implemented by `for n=1:numel(nu)`.
- Lines 33-34: Save results; implemented by `savefig(gcf,'xix_q_nutation_ensemble_b1_r.fig')`.

### Control flow inferred from the code

- Line 29: `for` loop over `n=1:numel(nu)`.

### Key state/data transformations

- Lines 14: computes `nu` using `nu=1e6*[6.8 9.6 13.5 17.5 25 36]`.
- Lines 17: computes `srt` using `srt=1e-3*[0.051 0.051 0.102 0.153 0.153 0.306]`.

### Local helper functions

- Line 38: `xix_field_profile_b1_r()` — `function xix_field_profile_b1_r(nu,srt)`. Q-band magnet
  - Representative operation: `sys.magnet=1.2142`.
  - Representative operation: `sys.isotopes={'E','1H'}`.

## Implementation structure

- Simulation of nutation frequency dependence of XiX DNP
- field profiles in the steady state with electron-proton
- distance and electron Rabi frequency ensembles.
- Calculation time: minutes
- Nutation frequencies, Hz
- Shot repetition times, seconds
- Get the figure started
- Plot the curves
- Save results
- Q-band magnet
- Electron and proton
- Zeeman interactions (g-tensor for trityl, ppm guess for 1H)

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `kfigure()`, `kxlabel()`, `kylabel()`, `kzlabel()`, `view()`, `xlim()`, `ylim()`, `zlim()`, `set()`, `xix_field_profile_b1_r()`, `srt()`, `savefig()`, `gaussleg()`, `r1n_dnp()`, `create()`, `basis()`.
