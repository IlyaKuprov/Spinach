# kernel/assume.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/kernel/assume.m`
- Signature: `spin_system=assume(spin_system,assumptions,retention)`
- Total lines: 769

## Purpose

Sets case-specific assumptions for various simulation contexts. This function determines the behaviour of the Hamiltonian generation func- tion and should be called before the Hamiltonian is requested. The function text is self-explanatory -interaction strength parameters are set in each section according to the physical requirements of of each specific simulation context. Syntax: spin_system=assume(spin_system,assum

## Physical / mathematical content

- This file belongs to the `kernel` part of Spinach. Its role should be read together with nearby files in the same directory, which usually share a common physical regime or infrastructure purpose.
- Quadrupolar physics is relevant: nuclei with spin > 1/2 interact with the electric field gradient tensor, introducing second-rank anisotropy, asymmetry, and overtone or MQ phenomena.
- The effective hardware model is a weakly anharmonic oscillator. Duffing nonlinearity breaks equal level spacing and allows qubit-like addressability within a truncated bosonic ladder.

## Numerical / algorithmic content

- The file contains an explicit `grumble(...)` validator, which is Spinach convention for front-loading dimension, type, and regime checks before expensive linear-algebra work begins.
- The file also defines local helper function(s): `grumble()`. This usually means the public entry point is supported by tightly coupled validation or helper logic kept private to the file.

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

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `grumble()`, `any()`, `ismember()`, `report()`, `strcmp()`, `specialised()`, `isfield()`, `elseif()`, `all()`, `cellfun()`, `exist()`, `ischar()`.
