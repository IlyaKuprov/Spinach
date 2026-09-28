# examples/optimal_control/state_transfer_s2m.m

- Signature: `state_transfer_s2m()`

## Purpose

Transfer coherence from a two-proton singlet state to a nearby carbon in a setting typical of parahydrogenation experiments. The example uses LBFGS-GRAPE as described at http://dx.doi.org/10.1016/j.jmr.2011.07.023. The source reports a terminal fidelity of 50% and a calculation time of minutes.

## Physical / mathematical content

- The system contains two `1H` and two `13C` spins at a magnetic field of 14.1, with scalar Zeeman values `{1.5, 2.0, 30.0, 40.0}` and specified scalar couplings of 7.0, 150, 150, and 50.
- The initial state is a singlet on spins 1 and 2; the target is `Lz` on spin 4. Each is normalised by its vector 2-norm.
- Proton and carbon `Lx` and `Ly` operators provide four controls. The drift Hamiltonian uses the `nmr` assumption and transmitter offsets of 1050 and 5285 for `1H` and `13C`, respectively.

## Numerical / algorithmic content

- The calculation uses the `sphten-liouv` formalism with no basis approximation. Its controls comprise 100 slices of duration `1.5e-4`, pulse-power levels `2*pi*[460 480 500 520 540]`, and `NS` and `SNS` penalties weighted 0.1 and 10.
- The control method is `lbfgs`, with a maximum of 100 iterations. Optimisation starts from a random `4`-by-`100` guess and calls `fmaxnewton` with `@grape_xy`.
- The resulting pulse is scaled by the mean power level. A test simulation applies it with `shaped_pulse_xy` using `expv-pwc` propagation and reports `real(rho_targ'*rho)`.

## Implementation structure

- Create the spin system and basis, construct and normalise the initial and target states, then assemble the control operators and offset-adjusted drift Hamiltonian.
- Configure the controls and optimisation plots, run `optimcon` and pulse optimisation, then scale the pulse and test it by simulation.
