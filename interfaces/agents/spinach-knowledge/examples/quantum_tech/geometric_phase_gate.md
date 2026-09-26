# examples/quantum_tech/geometric_phase_gate.m

- Signature: `geometric_phase_gate()`

## Purpose

Geometric phase gate between two trapped-ion qubits driven by a state-dependent optical dipole force, as demonstrated by Leibfried et al. (Nature 422, 412 (2003)) with two 9Be+ ions, simulated for a three-ion chain in the same trap: the global force beam pair detuned by delta from the stretch mode pushes the outer ions in opposite directions, while the middle ion has zero amplitude in the stretch mode eigenvector; the outer ions acquire a state-dependent geometric phase and the middle ion is a spectator that still couples to the far-detuned centre-of-mass and Egyptian modes. Calculation time: seconds

## Physical / mathematical content

- Ion qubits are spin-1/2 particles, the three axial normal modes of the chain (stretch at 6.1 MHz, centre-of-mass at stretch/sqrt(3), Egyptian at sqrt(29/5) times centre-of-mass, James, Appl. Phys. B 66, 181 (1998)) are bosonic modes `V6`, `V3`, `V3`; the spin-dependent force is a static `inter.modes.longitudinal` coupling in the frame rotating with the force drive, and the mode frequencies are laboratory values brought into that frame by `parameters.mode_offset` of the device context under the `spin-phonon` assumption set.
- The stretch mode force constant `kappa=delta/4` closes the phase-space loop at `T=1/delta` with a pi/2 differential geometric phase between the outer ions; the force constants of the other modes follow from the normal mode vectors and the 1/sqrt(frequency) scaling of the zero-point motion.
- Observables are normalised by the overlap with the unit state so that they are true expectation values; a second global pi/2 pulse turns the phase-gated state into a GHZ-type state of the outer ions, and their parity after an analysis pi/2 pulse of variable phase oscillates with full contrast, the measurement of Leibfried et al.

## Numerical / algorithmic content

- Basis set `IK-SBS` with `bas.inter_level=[2 3 2]`: boson-boson, spin-boson, and spin-spin coupling graphs traced separately; pure spin-spin correlations up to order two are kept inside the spin-boson subgraphs, which carries the entangling phase. The complete basis for this system has 186624 states.
- The pulse sequence is a local function passed to `device`: a global pi/2 pulse, an `evolution` trajectory over one closed loop of the stretch mode (`dt=gate_time/(npoints-1)`), a second global pi/2 pulse, and a loop of analysis pulses for the parity scan.
