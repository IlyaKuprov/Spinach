# experiments/hyperpol/dnp_freq_scan.m

- MATLAB implementation: [experiments/hyperpol/dnp_freq_scan.m](https://github.com/IlyaKuprov/Spinach/blob/main/experiments/hyperpol/dnp_freq_scan.m)

- Signature: `dnp=dnp_freq_scan(spin_system,parameters,H,R,K)`

## Purpose and physical scope

This routine computes steady-state DNP detection signals across a supplied microwave-frequency-offset vector. It does not sweep the static magnetic field or propagate a time-domain pulse sequence. Microwave drive, electron offset and electron-nuclear couplings act through the spin system and operators supplied by the caller; the function does not construct hyperfine tensors. It is not an ESEEM/ENDOR sequence or image-reconstruction routine.

## Inputs

- `parameters.mw_pwr`: scalar microwave power in radians per second.
- `parameters.mw_frq`: vector of frequency offsets in radians per second, specified relative to `parameters.g_ref`.
- `parameters.g_ref`: scalar reference g-factor (dimensionless).
- `parameters.rho0`: thermal-equilibrium state; `parameters.coil`: one detection-state vector or a horizontal stack.
- `parameters.mw_oper`: microwave irradiation operator. Liouville methods also require `parameters.ez_oper`, the electron `Lz` operator.
- `parameters.method`: `'lvn-backs'`, `'lvn-gmres'`, `'fp-backs'` or `'fp-gmres'`. The Fokker-Planck methods additionally require integer `parameters.nphases`, the microwave-phase grid size.
- H, R and K: Hamiltonian, relaxation and kinetics matrices supplied by the context function.

## Calculation and output axes

The Liouville-space paths convert offsets using the reference g-factor and magnet field, build the driven generator from H, R and K, then solve a steady-state linear system for each offset and project onto the `coil` states. The Fokker-Planck paths instead represent microwave phase on the `nphases` Fourier grid before solving. The output dnp has shape [numel(`parameters.mw_frq`), size(`parameters.coil`,2)]: rows follow the input frequency-vector order; columns follow the detection-state stack. The calculation returns expectation values, not a measured spectrum.

## Model limits

The source says R must not be thermalised for this calculation (`inter.equilibrium`='zero'). The supported formalisms are `sphten-liouv` and `zeeman-liouv`. The two solver choices in each method family are direct backslash and GMRES.

## Source-coded numerical example

`examples/dnp_sol/solid_effect_freq_scan_1.m` uses `parameters.mw_pwr`=2*pi*100e3 radians per second and two 100-point offset bands written as 2*pi*[linspace(144.0,145.5,100), linspace(14.0,15.5,100)]*1e6 radians per second. It sets a fixed orientation [pi/4 pi/5 pi/6] and uses `'lvn-backs'` for a single-crystal 15N-urea example. This is an example configuration, not a measurement or reported result.

## Source and attribution

- Source: `experiments/hyperpol/dnp_freq_scan.m`
- <https://spindynamics.org/wiki/index.php?title=dnp_freq_scan.m>
- Source attributions: ilya.kuprov@weizmann.ac.il; alexander.karabanov@nottingham.ac.uk; walter.kockenberger@nottingham.ac.uk; mariagrazia.concilio@sjtu.edu.cn
