# examples/optimal_control/state_transfer_coop.m

- Signature: `state_transfer_coop()`
- Source: [`examples/optimal_control/state_transfer_coop.m`](https://github.com/IlyaKuprov/Spinach/blob/main/examples/optimal_control/state_transfer_coop.m)

## Purpose

An optimal-control design for two pulses optimised cooperatively in a quadrupolar `14N` spin at a fixed orientation and power level. The stated design aim is for the combined outcomes to contain the target state without impurities; the source describes the calculation time as minutes. This is an objective, not a reported convergence result.

## Spin model and state transfer

The example uses a glycine nitrogen quadrupole interaction constructed with `eeqq2nqi(1.18e6,0.53,1,[1 2 3])`, a `14N` Zeeman scalar of 32.4, and `sys.magnet=14.1`. It uses the full spherical-tensor Liouville-space basis with no approximation. The initial and target states are the individually normalised `T1,0` and `T2,0` states of `14N`. The drift is formed from the isotropic and orientation-dependent quadrupolar Hamiltonian contributions, then transformed to the nitrogen carrier frame.

## Cooperative optimisation and diagnostics

The controls are `Lx` and `Ly`; the configured pulse has 100 slices of 10⁻⁷ s each (10 μs total), constant amplitude, and power level `2*pi*50e3`. The source supplies a random 2-by-100 initial guess, selects L-BFGS with a 100-iteration termination limit, and calls `fmaxnewton(spin_system,@grape_coop,guess)`. It requests coherence-order and phase-control plots. Finally, it prints the two final trajectory states, their average, and the target for comparison. The source contains no recorded numerical fidelity.
