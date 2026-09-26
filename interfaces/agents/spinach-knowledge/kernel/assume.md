# kernel/assume.m

- Signature: `spin_system=assume(spin_system,assumptions,retention)`

## Purpose

Sets case-specific assumptions for various simulation contexts. This function determines the behaviour of the Hamiltonian generation func- tion and should be called before the Hamiltonian is requested. The function text is self-explanatory -interaction strength parameters are set in each section according to the physical requirements of of each specific simulation context. Syntax: spin_system=assume(spin_system,assum

## Physical / mathematical content

- Quadrupolar physics is relevant: nuclei with spin > 1/2 interact with the electric field gradient tensor, introducing second-rank anisotropy, asymmetry, and overtone or MQ phenomena.
- The effective hardware model is a weakly anharmonic oscillator. Duffing nonlinearity breaks equal level spacing and allows qubit-like addressability within a truncated bosonic ladder.

## Numerical / algorithmic content

## Parameters / inputs

- assumptions -'nmr' for high-field NMR)
- 'esr' for electron rotating frame ESR
- 'deer' for DEER spectroscopy
- 'deer-zz' for DEER spectroscopy with electron
- flip-flop terms removed
- 'labframe' for full laboratory frame simulation
- with all Hamiltonian terms retained;
- bosonic modes are allowed and stay in
- the laboratory frame with all of their
- interaction terms retained
- 'qnmr' for quadrupolar NMR with numerical
- rotating frames: spin-1/2 particles
- will be in the rotating frame but
- spin>1/2 particles initially in the
- laboratory frame
- 'cavity' for cavity QED: spins and bosonic
- modes in a common rotating frame with
- the rotating wave approximation, mode
- energies to be built as detunings from
- the carrier frequency, exchange terms
- keeping flip-flop components only,
- anharmonicity, Kerr, and dispersive
- terms in full; longitudinal and modu-
- lation terms are disallowed because
- they average out
- 'spin-phonon' for spins in their usual rotating
- frames with bosonic modes in the la-
- boratory frame: electron and nuclear
- terms as in the 'esr' set, electron-
- mode exchange terms dropped as non-
- secular, mode-mode and nucleus-mode
- exchange terms retained in full,
- longitudinal, dispersive, and modula-
- tion terms and all diagonal mode
- terms retained
- retention -'zeeman' drops all spin-spin interactions
- 'couplings' drops all Zeeman interactions

## Outputs

- the function updates the spin_system object

## Implementation structure

- Sets case-specific assumptions for various simulation contexts. This
- function determines the behaviour of the Hamiltonian generation func-
- tion and should be called before the Hamiltonian is requested. The
- function text is self-explanatory -interaction strength parameters
- are set in each section according to the physical requirements of
- of each specific simulation context. Syntax:
- spin_system=assume(spin_system,assumptions,retention)
- assumptions -'nmr' for high-field NMR)
- 'esr' for electron rotating frame ESR
- 'deer' for DEER spectroscopy
- 'deer-zz' for DEER spectroscopy with electron
- flip-flop terms removed
