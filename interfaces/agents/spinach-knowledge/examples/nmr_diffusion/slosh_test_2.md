# examples/nmr_diffusion/slosh_test_2.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/examples/nmr_diffusion/slosh_test_2.m`
- Signature: `slosh_test_2()`
- Total lines: 37

## Purpose

Probability density sloshing around a harmonic oscillator in the presence of a gravitational pull twards the left. Calculation time: seconds.

## Physical / mathematical content

- Diffusion examples. The dominant mathematics is diffusion or advection-diffusion PDE propagation, sometimes with additional spin phase accumulation under gradients.

## Numerical / algorithmic content

- Time propagation is explicit. In Spinach this usually means repeated application of matrix exponentials or propagator factorizations to density operators or state vectors in Hilbert/Liouville/Fokker-Planck space.

## Code-derived implementation details

### Comment-guided execution stages

- Lines 10-11: Set oscillator parameters; implemented by `parameters.frc_cnst=2e3`.
- Lines 17-18: Get the Hamiltonian; implemented by `[H,~,xgrid]=oscillator(parameters)`.
- Lines 20-21: Get the initial state; implemented by `psi=exp(-50*(xgrid-0.1).^2)`.
- Lines 23-24: Get the propagator; implemented by `P=expm(full(-1i*H*0.001))`.
- Lines 26-27: Run the evolution; implemented by `kfigure()`.

### Control flow inferred from the code

- Line 28: `for` loop over `n=1:1000`.

### Key state/data transformations

- Lines 11: computes `parameters.frc_cnst` using `parameters.frc_cnst=2e3`.
- Lines 12: computes `parameters.par_mass` using `parameters.par_mass=1.00`.
- Lines 13: computes `parameters.box_size` using `parameters.box_size=2.00`.
- Lines 14: computes `parameters.n_points` using `parameters.n_points=100`.
- Lines 15: computes `parameters.grv_cnst` using `parameters.grv_cnst=250`.
- Lines 18: computes `[H,~,xgrid]` using `[H,~,xgrid]=oscillator(parameters)`.
- Lines 21: computes `psi` using `psi=exp(-50*(xgrid-0.1).^2)`.
- Lines 24: computes `P` using `P=expm(full(-1i*H*0.001))`.
- Lines 33: computes `drawnow; hold off; pause(0.001); psi` using `drawnow; hold off; pause(0.001); psi=P*psi`.

## Implementation structure

- Probability density sloshing around a harmonic oscillator
- in the presence of a gravitational pull twards the left.
- Calculation time: seconds.
- Set oscillator parameters
- Get the Hamiltonian
- Get the initial state
- Get the propagator
- Run the evolution

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `oscillator()`, `kfigure()`, `kxlabel()`, `kylabel()`, `pause()`.
