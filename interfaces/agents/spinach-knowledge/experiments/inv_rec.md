# experiments/inv_rec.m

[MATLAB source](https://github.com/IlyaKuprov/Spinach/blob/main/experiments/inv_rec.m) · [Spinach Wiki](https://spindynamics.org/wiki/index.php?title=inv_rec.m)

## Purpose

`inv_rec` is a liquid-NMR inversion-recovery sequence that generates a family of free-induction decays over a relaxation-delay trajectory. It starts from isotropic thermal equilibrium, applies a 180° inversion pulse, evolves the state under relaxation (and kinetics), applies a 90° pulse at each trajectory point, and detects an FID. The relaxation superoperator must be thermalised. This is a sequence implementation, not a measured relaxation result.

## Inputs and units

The signature takes `H`, `R`, and `K` and constructs `L=H+1i*R+1i*K`; unlike the imaging routines, it does not take `F` or gradient operators.

- `parameters.sweep`: spectral sweep width in Hz; acquisition step is `1/sweep` seconds.
- `parameters.npoints`: number of acquired FID points.
- `parameters.spins`: cell array of spin-name strings that selects the pulse and detection nucleus; source examples include `{'1H'}` and `{'13C'}`.
- `parameters.max_delay`: maximum relaxation-evolution duration, in seconds.
- `parameters.n_delays`: positive integer number of relaxation-evolution steps spanning `max_delay` (step interval `max_delay/n_delays`).

The pulse/detection states are formed from `L+` of the selected first spin: `Ly` generates the 180° and 90° rotations, while `coil_state(spin_system,'L+',parameters.spins{1},'exact')` provides the detection observable. There is no explicit coherence-order filter.

## Detection and return

`fids` contains the acquired FIDs as columns, with one column per relaxation-trajectory state; the row dimension is the FID time-point dimension (`npoints`). The trajectory starts at zero delay, and each propagation interval is `max_delay/n_delays`. The function returns the matrix only, not explicit time or delay-axis vectors. No DOI or numerical experimental data are supplied in the source page.
