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

## Code-derived implementation details

### Comment-guided execution stages

- Lines 69-70: Check consistency; implemented by `if nargin==2`.
- Lines 78-79: Store the assumption specification; implemented by `spin_system.inter.assumptions=assumptions`.
- Lines 81-82: Preallocate the arrays; implemented by `spin_system.inter.zeeman.strength=cell(spin_system.comp.nspins,1)`.
- Lines 86-87: Set the approximations; implemented by `switch assumptions`.
- Lines 91-92: Disallow particles that are not spins; implemented by `if any(ismember({'C','V','T'},spin_system.comp.types),'all')`.
- Lines 96-97: Do the reporting; implemented by `report(spin_system,'generic high-field NMR assumption set:')`.
- Lines 107-108: Process Zeeman interactions; implemented by `for n=1:spin_system.comp.nspins`.
- Lines 110-111: All Zeeman interactions should be secular; implemented by `spin_system.inter.zeeman.strength{n}='secular'`.
- Lines 115-116: Process couplings; implemented by `for n=1:spin_system.comp.nspins`.
- Lines 118-119: Giant spin terms should be secular; implemented by `spin_system.inter.giant.strength{n}='secular'`.
- Lines 121-122: For spin-spin couplings, it depends; implemented by `for k=1:spin_system.comp.nspins`.
- Lines 125-126: Couplings between spins of the same type should be secular; implemented by `spin_system.inter.coupling.strength{n,k}='secular'`.
- Lines 130-131: Couplings between spins of different types should be weak; implemented by `spin_system.inter.coupling.strength{n,k}='zz'`.
- Lines 145-146: Do the reporting; implemented by `report(spin_system,'generic high-field EPR/DEER assumption set:')`.
- Lines 162-163: Electron Zeeman interactions should be secular; implemented by `spin_system.inter.zeeman.strength{n}='secular'`.
- Lines 167-168: Full Zeeman tensors should be used for nuclei; implemented by `spin_system.inter.zeeman.strength{n}='full'`.
- Lines 183-184: Couplings from electrons to nuclei should be left-secular; implemented by `spin_system.inter.coupling.strength{n,k}='z*'`.
- Lines 188-189: Couplings from nuclei to electrons should be right-secular; implemented by `spin_system.inter.coupling.strength{n,k}='*z'`.

### Control flow inferred from the code

- Line 70: conditional branch on `nargin==2`.
- Line 87: dispatches on `assumptions`; cases `{'nmr'}`, `{'esr','deer'}`, `{'deer-zz'}`, `{'labframe'}`.
- Line 92: conditional branch on `any(ismember({'C','V','T'},spin_system.comp.types),'all')`.
- Line 108: `for` loop over `n=1:spin_system.comp.nspins`.
- Line 116: `for` loop over `n=1:spin_system.comp.nspins`.
- Line 122: `for` loop over `k=1:spin_system.comp.nspins`.
- Line 123: conditional branch on `strcmp(spin_system.comp.isotopes{n},spin_system.comp.isotopes{k})`.
- Line 141: conditional branch on `any(ismember({'C','V','T'},spin_system.comp.types),'all')`.
- Line 159: `for` loop over `n=1:spin_system.comp.nspins`.
- Line 160: conditional branch on `strcmp(spin_system.comp.isotopes{n}(1),'E')`.
- Line 174: `for` loop over `n=1:spin_system.comp.nspins`.
- Line 180: `for` loop over `k=1:spin_system.comp.nspins`.
- Line 181: conditional branch on `strcmp(spin_system.comp.isotopes{n}(1),'E')&&(~strcmp(spin_system.comp.isotopes{k}(1),'E'))`.
- Line 208: conditional branch on `any(ismember({'C','V','T'},spin_system.comp.types),'all')`.

### Key state/data transformations

- Lines 79: computes `spin_system.inter.assumptions` using `spin_system.inter.assumptions=assumptions`.
- Lines 82: computes `spin_system.inter.zeeman.strength` using `spin_system.inter.zeeman.strength=cell(spin_system.comp.nspins,1)`.
- Lines 83: computes `spin_system.inter.giant.strength` using `spin_system.inter.giant.strength=cell(spin_system.comp.nspins,1)`.
- Lines 84: computes `spin_system.inter.coupling.strength` using `spin_system.inter.coupling.strength=cell(spin_system.comp.nspins)`.
- Lines 111: computes `spin_system.inter.zeeman.strength{n}` using `spin_system.inter.zeeman.strength{n}='secular'`.
- Lines 119: computes `spin_system.inter.giant.strength{n}` using `spin_system.inter.giant.strength{n}='secular'`.
- Lines 126: computes `spin_system.inter.coupling.strength{n,k}` using `spin_system.inter.coupling.strength{n,k}='secular'`.
- Lines 323: computes `spin_system.inter.modes.strength.frqs` using `spin_system.inter.modes.strength.frqs='full'`.
- Lines 324: computes `spin_system.inter.modes.strength.anharms` using `spin_system.inter.modes.strength.anharms='full'`.
- Lines 325: computes `spin_system.inter.modes.strength.exchange` using `spin_system.inter.modes.strength.exchange='strong'`.
- Lines 326: computes `spin_system.inter.modes.strength.kerr` using `spin_system.inter.modes.strength.kerr='full'`.
- Lines 327: computes `spin_system.inter.modes.strength.longitudinal` using `spin_system.inter.modes.strength.longitudinal='full'`.
- Lines 328: computes `spin_system.inter.modes.strength.dispersive` using `spin_system.inter.modes.strength.dispersive='full'`.
- Lines 329: computes `spin_system.inter.modes.strength.coupling_mod` using `spin_system.inter.modes.strength.coupling_mod='full'`.
- Lines 330: computes `spin_system.inter.modes.strength.zeeman_mod` using `spin_system.inter.modes.strength.zeeman_mod='full'`.
- Lines 519: computes `report(spin_system,' rotating frame approximation for S` using `report(spin_system,' rotating frame approximation for S=1/2 nuclei,')`.
- Lines 522: computes `report(spin_system,' secular Zeeman terms for S` using `report(spin_system,' secular Zeeman terms for S=1/2 nuclei,')`.
- Lines 524: computes `report(spin_system,' secular terms for couplings between S` using `report(spin_system,' secular terms for couplings between S=1/2 nuclei,')`.

### Local helper functions

- Line 755: `grumble()` — `function grumble(assumptions,retention)`. Humanity has advanced, when it has advanced, not because it has been sober, responsible, and cautious, but because it has been
  - Representative operation: `if ~ischar(assumptions)`.
  - Representative operation: `error('assumptions argument muct be a character string.')`.

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
