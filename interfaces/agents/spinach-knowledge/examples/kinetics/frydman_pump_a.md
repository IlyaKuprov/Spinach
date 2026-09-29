# examples/kinetics/frydman_pump_a.m

- MATLAB implementation: [examples/kinetics/frydman_pump_a.m](https://github.com/IlyaKuprov/Spinach/blob/main/examples/kinetics/frydman_pump_a.m)

- Callable as the no-argument MATLAB function `frydman_pump_a()`. It creates the spin system and calls `liquid` with the local `frydman_pump` sequence callback.

## Purpose and model

The source identifies this as Lucio Frydman's water-exchange spin-lock pump, Figure 2, and cites https://doi.org/10.1016/j.jmr.2021.107083. The model has a `1H`, `15N`, `13C`, `13C` peptide core plus 20 added `1H` water spins (24 spins total). `sys.magnet` is set to 11.7; the source labels this a magnet field but gives no unit. All listed chemical shifts are zero. The explicitly assigned H-N scalar coupling is -45 Hz.

The `t1_t2` relaxation model uses diagonal retention. Source-estimated T1 and T2 values are given in seconds: H 0.2722, N 0.8, C-alpha 2, carbonyl C 2, and water 0.2994. The source sets `inter.temperature=298` without stating a unit and uses an isotropic equilibrium model (`IME`). The basis uses `sphten-liouv`, `IK-1`, scalar-coupling connectivity, `inter_level=4`, and `prox_level=1`. Intermolecular exchange entries are source-labelled in Hz: 10 between the amide H and the first water proton, and 1e4 in the water-proton block.

## Sequence and observables

The NMR call selects `1H`, `15N`, and `13C`, with a 100 ms CP duration and 100 points. In the callback, the source removes Zeeman terms for these nuclei and treats H-N, N-C, and H-C couplings as strong when constructing the effective spin-lock Hamiltonian. Starting from `equilibrium(spin_system)`, it applies a `+pi/2` N-y step, destroys the remaining components with `homospoil`, then applies `+pi/2` H-y and N-y steps. The CP trajectory is propagated with `krylov` using a 1 ms time step and `cp_npt-1` steps.

The reported observables are real projections of the trajectory onto H-z, H-x, N-z, and N-x operators for the peptide-bond H and N. Two subplots show H Z/X and N Z/X expectation values; the axis is labelled time in ms. The source header estimates seconds of calculation time, not a measured runtime.

## Scope

The script supplies a simulation setup and identifies the cited figure, but no trajectory values or comparison with Figure 2. Accordingly, it does not support claims about quantitative reproduction, transfer efficiency, or a measured exchange result.