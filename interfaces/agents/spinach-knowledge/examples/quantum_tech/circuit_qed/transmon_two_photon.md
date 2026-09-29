# examples/quantum_tech/circuit_qed/transmon_two_photon.m

- Source: [transmon_two_photon.m](https://github.com/IlyaKuprov/Spinach/blob/main/examples/quantum_tech/circuit_qed/transmon_two_photon.m)
- Signature: `transmon_two_photon()`

## Purpose

Find a shaped pulse that drives a four-level Duffing transmon from its ground state to its second excited state by a two-photon transition, while the carrier is off resonance from both adjacent one-photon transitions. This is a simulated control trajectory, not a measured device transfer.

## Hamiltonian and rotating frame

The system declares the four-state transmon mode as `T4`, with detuning `inter.modes.frqs={100e6}` and anharmonicity `inter.modes.anharms={-200e6}` (frequency settings in Hz: 100 MHz and −200 MHz). In the Duffing/Fock description, the mode drift has the form `H ∝ Δ a†a + (α/2) a†a†aa`, truncated to four levels; here Δ=100 MHz and α=−200 MHz. Its 0↔1 and 1↔2 transition offsets are therefore +100 MHz and −100 MHz, so the carrier sits halfway between them. There is no separately parameterised resonator mode or transmon–resonator coupling in this script.

The example constructs `H=hamiltonian(assume(spin_system,'cavity'))` using the `zeeman-hilb` basis and no basis approximation. In this Spinach cavity-QED frame, mode energies are carrier detunings and the rotating-wave approximation is used; the example's Duffing anharmonicity is retained. The oscillator quadratures are built from the transmon ladder operators as `Cx=(C+A)/2` and `Cy=i(C−A)/2`.

## Control sequence and objective

The source and target are the normalised `BL1` and `BL3` states (ground and second excited levels in the source comments). GRAPE uses the fixed drift `H` and the two quadrature controls `Cx`, `Cy`; `fmaxnewton` optimises the pulse with `grape_xy` from a smooth 200-sample initial guess. The target is state transfer to the second excited level via the virtual intermediate state, not a resonant one-photon drive.

The optimised pulse is then propagated slice by slice. The plotted observable is the population in each of the four levels versus accumulated pulse time, labelled in nanoseconds. It is the model's simulated trajectory; this script does not calculate a device readout signal or report experimental hyperpolarisation/measurement data.

## Provenance

The source comments identify the model and parameters as Example 2 of the GRAPE_SCQ package and estimate minutes of calculation time. That is source-provided context, not a runtime measured for this note.
