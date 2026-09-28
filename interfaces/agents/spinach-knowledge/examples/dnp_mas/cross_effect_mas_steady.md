# examples/dnp_mas/cross_effect_mas_steady.m

- Signature: `cross_effect_mas_steady()`

## Purpose

Finds the periodic steady state of a single-crystal MAS DNP model, then analyses its trajectory over one rotor period and reports the proton enhancement. The example follows [Mentink-Vigier et al.](http://dx.doi.org/10.1016/j.jmr.2015.07.001), noting that Spinach uses different rotation conventions; the stated runtime is minutes.

## Physical and numerical setup

The system comprises two electrons and one proton at 9.394 T, with anisotropic electron g tensors, electron-electron and electron-proton couplings, and Nottingham relaxation. The specified electron/nuclear relaxation rates, 100 K temperature, Di Bari equilibrium, and secular relaxation retention are used with the complete spherical-tensor Liouville basis. A magnetic-frame ESR rotor stack is generated for the selected electron and proton spins with the given rotor axis and crystal orientation, using rank 3000.

For each rotor interval, the code combines the stack Hamiltonian with electron microwave-drive and offset terms and the relaxation superoperator. It composes the interval propagators over one rotor period and obtains the periodic steady state with `steady(...,'newton')`. Starting from that state, it steps through one period, plots level populations with `trajan`, and computes the proton enhancement as the time-mean proton Lz signal divided by the thermal-equilibrium proton Lz signal.
