# examples/optimal_control/state_transfer_m2s.m

- Signature: `state_transfer_m2s()`

## Purpose

A transfer of coherence from longitudinal magnetization into a two-spin singlet state in allyl pyruvate with a distribution of B1 powers and transmitter offsets. XX and YY components of the singlet dephase rapidly in this system, and are therefore dropped from the target state specification. Calculation time: many hours.

## Physical / mathematical content

- The initial state is longitudinal magnetization on all spins; the target is the negative two-spin `Lz`–`Lz` state on `Hb` and `Hc`.
- Proton `Lx` and `Ly` operators provide the controls. The drift Hamiltonian includes a 2670 Hz transmitter offset, with an additional five-point offset distribution from −10 to 10 Hz.

## Numerical / algorithmic content

- GRAPE optimisation uses L-BFGS with five pulse-power levels, `2*pi*[480 490 500 510 520]`, 300 slices of 1 ms, `NS` and `SNS` penalties weighted `[1 10]`, and a maximum of 3000 iterations.
- A test simulation applies the optimised pulse, destroys transverse coherence with homospoiling, and reports `Re[<target|rho(T)>]`. Pulse-acquire spectra before and after the pulse are exponentially apodised and Fourier transformed.

## Implementation structure

- Load the allyl pyruvate proton spin system, remove the methyl group, set the magnetic field to 11.7464 T (500.13 MHz), and build an unrestricted spherical-tensor Liouville-space basis.
- Construct the initial and target states, drift and control operators, ensemble, and initial pulse guess; then optimise with `fmaxnewton(spin_system,@grape_xy,pulse)`.
- Simulate and plot the initial pulse-acquire spectrum alongside the spectrum obtained after applying the optimised pulse.
