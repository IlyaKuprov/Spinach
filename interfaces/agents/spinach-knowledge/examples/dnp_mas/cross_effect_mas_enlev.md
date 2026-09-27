# examples/dnp_mas/cross_effect_mas_enlev.m

- Signature: `cross_effect_mas_enlev()`

## Purpose

Plots the spin energy levels over one MAS rotor period for a single-crystal model associated with the cross-effect DNP example of [Mentink-Vigier et al.](http://dx.doi.org/10.1016/j.jmr.2015.07.001). The source notes that Spinach uses different rotation conventions from the paper; its stated runtime is milliseconds.

## Physical and numerical setup

The system has two electron spins and one proton at 9.394 T, with anisotropic electron g tensors and electron-electron and electron-proton couplings. It uses the complete Zeeman-Hilbert basis. A lab-frame rotor stack is generated at 12.5 kHz for the selected electron and proton spins, with the specified rotor axis, crystal orientation, and magnetic-frame MAS convention.

For each rotor-stack Hamiltonian, the code computes and sorts its real eigenvalues, converts them to GHz, and plots the lowest eight levels against time over one rotor period in three panels. The stack diagonalizations are run with `parfor`. This example is an energy-level analysis; it does not propagate a density operator or include relaxation or microwave-drive terms.
