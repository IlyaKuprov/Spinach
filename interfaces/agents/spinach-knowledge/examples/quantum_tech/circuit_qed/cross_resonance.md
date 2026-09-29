# examples/quantum_tech/circuit_qed/cross_resonance.m

- Signature: `cross_resonance()`
- Source: [cross_resonance.m](https://github.com/IlyaKuprov/Spinach/blob/main/examples/quantum_tech/circuit_qed/cross_resonance.m)

## Purpose and mechanism

The example illustrates the cross-resonance mechanism in two fixed-frequency transmons: drive the control transmon near the target’s transition frequency, and the static exchange coupling makes the target’s response depend on the control state. The two conditional target rotations differ, giving the ZX component used in cross-resonance gates. The source also notes that the unconditional response is larger at the selected setting, so an echo would be needed to cancel it; this file simulates a single drive interval, not that echo sequence.

## Model and drive

The model contains two three-level Duffing transmons (`T3`, `T3`). Their input frequencies are 5.5 and 6.0 GHz, their anharmonicities are −240 and −200 MHz, and their exchange coupling is 25 MHz. These are supplied through Spinach’s mode-frequency and coupling fields in Hz, which `create` converts to angular-frequency units internally. The laboratory-frame drift contains both Duffing-mode energies (their transition frequencies and self-anharmonicities) and the static mode-exchange interaction; the two drive operators act separately on the transmons.

The common local oscillator is `2π × 6.0002 GHz`, 200 kHz above the target’s declared bare frequency. The two channel amplitudes are 5.969026041820607×10^8 and 2.883982055995430×10^7 rad/s, with phases 0 and −0.39720756 rad. A flat-top Gaussian envelope has rise and fall centres at 30 and 120 ns and a 15 ns ramp. Thus the propagated drive interval is 150 ns. The source describes these calibrated drive settings as fixed reference values, not parameters derived or optimised in this script.

## Simulated sequence and observable

For each run, the target begins in its ground state while the control is prepared either in its ground or first excited state. The script applies the same 150 ns driven evolution to each initial ket, stepping the laboratory-frame Hamiltonian at 10 ps intervals. Every 75 steps (0.75 ns), it rotates the kets into the drive frame and records the target’s three Pauli expectation values. The output is the pair of conditional target Bloch trajectories, not a process-fidelity estimate or an experimental readout.

The checks in the source require norm conservation within 10^-6 and impose signed final-Bloch-component bounds for the control-ground and control-excited cases. They are script assertions, not results newly verified here. No relaxation channel, device measurement, or echoed-gate sequence is included.

The source header attributes its model and parameters to the matching example in the paraqeet package.
