# examples/quantum_tech/geometric_phase_gate.m

[Example source](https://github.com/IlyaKuprov/Spinach/blob/main/examples/quantum_tech/geometric_phase_gate.m)

- Signature: `geometric_phase_gate()`

## Physical model

The script simulates three spin-1/2 ion qubits coupled to the three axial phonon modes of a linear chain. Its comments place the gate in the context of the two ⁹Be⁺-ion experiment of Leibfried et al., *Nature* 422, 412 (2003); this code instead models three ions in the same trap. The centre ion has zero amplitude in the stretch-mode eigenvector, while the two outer ions are driven with opposite stretch-mode force signs. The centre ion is a spectator for that mode but remains coupled to the off-resonant centre-of-mass and Egyptian modes.

In the rotating frame of the force drive, the spin-dependent force is a static longitudinal spin-boson coupling, entered through `inter.modes.longitudinal` with Spinach's `Lz(a+a')/sqrt(2)` convention. The basis is `sphten-liouv` with the `IK-SBS` approximation and `inter_level=[2 3 2]`; it retains the spin/phonon correlations needed by this model. The device context is called with the `spin-phonon` assumption set. With three two-level spins and the `V6`, `V3`, `V3` mode cutoffs, the full Hilbert space has 432 states and the corresponding Liouville space has 186,624 operators. No decoherence or motional-loss parameters are specified, so this is not a hardware-fidelity calculation.

## Drive and simulated observables

The stretch mode is set to 6.1 MHz. The force is detuned from it by 26 kHz, so the drive frequency is 6.126 MHz; the gate interval is `1/delta`, approximately 38.46 µs, and the stretch-mode force constant is `delta/4`, or 6.5 kHz. The other axial frequencies follow the source formulas: centre-of-mass at stretch divided by `sqrt(3)`, and Egyptian at `sqrt(29/5)` times centre-of-mass. The trajectory has 1,561 sample points.

The first plot shows the three ion σx expectation values and stretch-mode occupation ⟨a†a⟩ during the simulated gate interval. The second is a calculated parity scan of the outer ions versus analysis-pulse phase after the gate. The script includes checks for the spectator coherence, return of the stretch mode to its ground state, and parity contrast; listing those checks does not establish that a run passed them. These are simulated observables, not measurements or a reported gate fidelity.

The source cites Leibfried et al. (2003) for the geometric gate and James, *Appl. Phys. B* 66, 181 (1998) for the three-ion mode frequencies.
