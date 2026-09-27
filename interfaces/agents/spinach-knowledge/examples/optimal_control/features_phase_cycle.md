# examples/optimal_control/features_phase_cycle.m

- Signature: `features_phase_cycle()`

## Purpose

Optimise a pulse for state-to-state transfer across scalar couplings in a hydrofluorocarbon fragment spin system. The initial state is ¹H Z-magnetisation; the target is ¹⁹F transverse magnetisation. A two-step phase cycle tests whether flipping the fluorine-channel phase also flips the resulting ¹⁹F magnetisation. The source notes a calculation time of minutes.

## Spin system and states

- Magnetic field: `sys.magnet=9.4`.
- Isotopes: `{'1H','13C','19F'}`; chemical shifts: `{0.0,0.0,0.0}` ppm.
- Scalar couplings: `inter.coupling.scalar{1,2}=140` Hz and `inter.coupling.scalar{2,3}=-160` Hz.
- Basis: `sphten-liouv`, with approximation `none`.
- Normalised initial state: `state(spin_system,{'Lz'},{1})`.
- Normalised target state: `state(spin_system,{'L+'},{3})+state(spin_system,{'L-'},{3})`.

## Control and optimisation

- The drift is `hamiltonian(assume(spin_system,'nmr'))`. The six control operators are `Lx` and `Ly` for each of ¹H, ¹³C and ¹⁹F, with channel map `[1;1;2;2;3;3]`.
- Pulse-power levels are `2*pi*linspace(0.8e3,1.2e3,5)` rad/s; slice durations are `2e-4*ones(1,50)`. The penalty is `SNS` with weight `100`.
- The optimisation method is `lbfgs`, with `max_iter=100`. The phase-cycle matrix is `[0 0 0 0 0; 0 0 0 pi pi]`; the source notes that only `0` and `pi` are supported.
- After `optimcon`, the initial guess is `randn(6,50)/3`. Optimisation calls `fmaxnewton(spin_system,@grape_xy,guess)`, and the resulting pulse is multiplied by `mean(control.pwr_levels)`.

## Phase-cycle tests

For each row of `control.phase_cycle`, the script takes the phase from column 4, rotates fluorine control rows `5:6` by that phase, and simulates the phased pulse with `shaped_pulse_xy(...,'expv-pwc')`. It then reports state information using `stateinfo(spin_system,rho,10)`.