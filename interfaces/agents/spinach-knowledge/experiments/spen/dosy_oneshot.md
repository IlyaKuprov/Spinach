# experiments/spen/dosy_oneshot.m

- Signature: `fid=dosy_oneshot(spin_system,parameters,H,R,K,G,F)`

## Purpose

Implements a one-shot DOSY pulse sequence.

## Physical / mathematical content

The routine combines spin evolution with the supplied Fokker–Planck gradient and diffusion/flow superoperators. It selects coherence orders during the pulse sequence and uses bipolar gradient intervals around the diffusion period.

## Numerical / algorithmic content

Starting from `parameters.rho0`, it applies 90-degree and 180-degree pulses about the Y operator for the first specified spin, alternates coherence-order selection with gradient and stabilization-delay evolution, evolves for the remaining diffusion interval, and acquires the FID on `parameters.coil`. The acquired signal uses dwell time `1/parameters.sweep` and `parameters.npoints-1` intervals.

## Parameters / inputs

- `parameters.rho0` — initial state
- `parameters.coil` — detection state
- `parameters.spins` — nuclei on which the sequence runs
- `parameters.g_amp` — gradient amplitude for diffusion encoding, T/m
- `parameters.g_dur` — pulse width of the gradient for diffusion encoding, s
- `parameters.kappa` — unbalancing factor for the bipolar gradients, with ratio `(1+kappa):(1-kappa)`
- `parameters.g_stab_del` — gradient stabilization delay, s
- `parameters.del` — diffusion delay, s
- `parameters.dims` — sample size, m
- `parameters.npts` — number of discretization points in the grid
- `parameters.diff` or `parameters.dxx` — spatially uniform or voxel-wise diffusion coefficients, respectively
- `parameters.npoints` — number of points in the acquired signal
- `parameters.sweep` — acquisition sweep width, Hz
- `H`, `R`, `K`, `G`, `F` — Fokker–Planck Hamiltonian, relaxation, kinetics, gradient, and diffusion/flow superoperators, respectively

## Outputs

- `fid` — free induction decay

## Reference

- [Spinach documentation](https://spindynamics.org/wiki/index.php?title=dosy_oneshot.m)
