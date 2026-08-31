# examples/quantum_tech/transmon_duffing_ladder.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/examples/quantum_tech/transmon_duffing_ladder.m`
- Signature: `transmon_duffing_ladder()`
- Total lines: 62

## Purpose

Duffing-model energy ladder of a weakly anharmonic transmon, showing how the transition frequencies separate as anharmo- nicity increases. Inspired by the transmon model of Koch et al., Phys. Rev. A 76, 042319 (2007). Calculation time: seconds

## Physical / mathematical content

- Quantum-technology examples. The files in this area model cavity QED, transmon qubits, NV centres, and related effective Hamiltonians. The recurring mathematics is finite-dimensional quantum dynamics with ladder operators, rotating-wave-style couplings, anharmonic oscillator terms, avoided crossings, and coherent control in coupled few-mode systems.
- The effective hardware model is a weakly anharmonic oscillator. Duffing nonlinearity breaks equal level spacing and allows qubit-like addressability within a truncated bosonic ladder.

## Numerical / algorithmic content

- An eigenvalue problem is solved or analysed, so the file is extracting spectra, stationary states, avoided crossings, or modal structure from the effective Hamiltonian or superoperator.

## Code-derived implementation details

### Comment-guided execution stages

- Lines 12-13: Magnet field; implemented by `sys.magnet=0`.
- Lines 15-16: Particle specification; implemented by `sys.isotopes={'T5'}`.
- Lines 18-19: Transmon mode frequency; implemented by `inter.modes.frqs={5.0e9}`.
- Lines 21-22: Formalism and basis; implemented by `bas.formalism='zeeman-hilb'`.
- Lines 25-26: Spinach housekeeping; implemented by `spin_system=create(sys,inter)`.
- Lines 29-30: Harmonic part from the declared frequency; implemented by `spin_system=assume(spin_system,'labframe')`.
- Lines 33-34: Anharmonicity operator; implemented by `K=operator(spin_system,'CCAA',1)`.
- Lines 36-37: Anharmonicity sweep range; implemented by `anharm=2*pi*linspace(-400e6,-50e6,80)`.
- Lines 39-40: Preallocate transition frequency array; implemented by `trans_frq=zeros(numel(anharm),4)`.
- Lines 42-43: Sweep the Duffing anharmonicity; implemented by `for n=1:numel(anharm)`.
- Lines 45-46: Build the Duffing Hamiltonian; implemented by `H=H0+(anharm(n)/2)*K`.
- Lines 48-49: Extract adjacent transition frequencies; implemented by `energies=sort(real(eig(full(H))))`.
- Lines 54-55: Plot the transition ladder; implemented by `kfigure(); plot(-anharm/(2*pi*1e6),trans_frq/1e9,'LineWidth',1.5)`.

### Control flow inferred from the code

- Line 43: `for` loop over `n=1:numel(anharm)`.

### Key state/data transformations

- Lines 13: computes `sys.magnet` using `sys.magnet=0`.
- Lines 16: computes `sys.isotopes` using `sys.isotopes={'T5'}`.
- Lines 19: computes `inter.modes.frqs` using `inter.modes.frqs={5.0e9}`.
- Lines 22: computes `bas.formalism` using `bas.formalism='zeeman-hilb'`.
- Lines 23: computes `bas.approximation` using `bas.approximation='none'`.
- Lines 26: computes `spin_system` using `spin_system=create(sys,inter)`.
- Lines 31: computes `H0` using `H0=hamiltonian(spin_system)`.
- Lines 34: computes `K` using `K=operator(spin_system,'CCAA',1)`.
- Lines 37: computes `anharm` using `anharm=2*pi*linspace(-400e6,-50e6,80)`.
- Lines 40: computes `trans_frq` using `trans_frq=zeros(numel(anharm),4)`.
- Lines 46: computes `H` using `H=H0+(anharm(n)/2)*K`.
- Lines 49: computes `energies` using `energies=sort(real(eig(full(H))))`.
- Lines 50: computes `trans_frq(n,:)` using `trans_frq(n,:)=diff(energies)/(2*pi)`.

## Implementation structure

- Duffing-model energy ladder of a weakly anharmonic transmon,
- showing how the transition frequencies separate as anharmo-
- nicity increases. Inspired by the transmon model of Koch et
- al., Phys. Rev. A 76, 042319 (2007).
- Calculation time: seconds
- Magnet field
- Particle specification
- Transmon mode frequency
- Formalism and basis
- Spinach housekeeping
- Harmonic part from the declared frequency
- Anharmonicity operator

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `create()`, `basis()`, `assume()`, `hamiltonian()`, `operator()`, `anharm()`, `trans_frq()`, `diff()`, `kfigure()`, `kxlabel()`, `kylabel()`, `ktitle()`, `klegend()`.
