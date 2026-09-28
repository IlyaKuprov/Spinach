# examples/optimal_control/state_transfer_coop.m

- Signature: `state_transfer_coop()`

## Purpose

Optimise two pulses cooperatively for state-to-state transfer in a quadrupolar 14N spin at a fixed orientation and power level, so that their combined outcomes contain the target state without impurities. Calculation time: minutes.

## Physical / mathematical content

- The initial state is the normalised 14N `T1,0` state; the target is the normalised `T2,0` state.
- The drift Hamiltonian is assembled at orientation `[1 2 3]` and transformed to a rotating frame. The control operators are `Lx` and `Ly`.
- The optimisation uses the cooperative GRAPE objective `grape_coop` with the limited-memory quasi-Newton method `lbfgs`.

## Numerical / algorithmic content

- The spin system has a 14.1 T field, a glycine 14N NQI coupling specified by `eeqq2nqi(1.18e6,0.53,1,[1.0 2.0 3.0])`, and a 14N chemical shift of 32.4. It uses the `sphten-liouv` formalism without basis approximation.
- The controls use a power level of `2*pi*50e3`, 100 slices of duration `10e-8`, an initial amplitude profile of ones, a random `2`-by-`100` initial guess, and a limit of 100 optimisation iterations.

## Implementation structure

- Create the spin system and basis, prepare and normalise the initial and target states, and construct the drift and control operators.
- Configure the cooperative optimisation, run `fmaxnewton` with `@grape_coop`, then print both final outcomes, their average, and the target state.
