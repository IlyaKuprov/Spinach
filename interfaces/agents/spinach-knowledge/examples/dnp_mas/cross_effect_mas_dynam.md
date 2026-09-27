# examples/dnp_mas/cross_effect_mas_dynam.m

- Signature: `cross_effect_mas_dynam()`

## Purpose

Tracks the spin-system trajectory through one MAS rotor period for a single crystal in a microwave-driven DNP simulation based on [Mentink-Vigier et al.](http://dx.doi.org/10.1016/j.jmr.2015.07.001). The example notes that Spinach uses different rotation conventions from the paper; the calculation is described as taking seconds.

## Physical and numerical setup

The model contains two electron spins and one proton at a 9.394 T field. It assigns anisotropic electron g tensors, an electron-electron coupling and an electron-proton coupling, and uses the Nottingham relaxation model with the specified electron/nuclear relaxation rates, 100 K temperature, Di Bari equilibrium, and secular relaxation retention. The spin system is represented in the complete spherical-tensor Liouville basis (no basis approximation).

A single-crystal rotor stack is generated in the magnetic frame with the listed rotor axis and orientation, electron/proton spin selection, and ESR settings. Starting from thermal equilibrium, the trajectory is stepped through the rotor-stack Hamiltonians at 12.5 kHz with a microwave drive amplitude parameter of 0.85 MHz and an electron offset of -400 MHz. The resulting density-operator trajectory is analysed with `trajan` at the level-population level.
