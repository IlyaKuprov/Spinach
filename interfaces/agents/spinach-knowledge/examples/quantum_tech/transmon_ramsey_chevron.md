# examples/quantum_tech/transmon_ramsey_chevron.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/examples/quantum_tech/transmon_ramsey_chevron.m`
- Signature: `transmon_ramsey_chevron()`
- Total lines: 86

## Purpose

Ramsey chevron of a three-level transmon in the Duffing ap- proximation. A nominal pi/2 pulse prepares a coherence, and detuning during free evolution produces Ramsey fringes. Calculation time: seconds

## Physical / mathematical content

- Quantum-technology examples. The files in this area model cavity QED, transmon qubits, NV centres, and related effective Hamiltonians. The recurring mathematics is finite-dimensional quantum dynamics with ladder operators, rotating-wave-style couplings, anharmonic oscillator terms, avoided crossings, and coherent control in coupled few-mode systems.
- The effective hardware model is a weakly anharmonic oscillator. Duffing nonlinearity breaks equal level spacing and allows qubit-like addressability within a truncated bosonic ladder.

## Numerical / algorithmic content

- Time propagation is explicit. In Spinach this usually means repeated application of matrix exponentials or propagator factorizations to density operators or state vectors in Hilbert/Liouville/Fokker-Planck space.

## Code-derived implementation details

### Comment-guided execution stages

- Lines 11-12: Magnet field; implemented by `sys.magnet=0`.
- Lines 14-15: Particle specification; implemented by `sys.isotopes={'T3'}`.
- Lines 17-18: Transmon in the rotating frame; implemented by `inter.modes.frqs={0}`.
- Lines 21-22: Formalism and basis; implemented by `bas.formalism='zeeman-hilb'`.
- Lines 25-26: Spinach housekeeping; implemented by `spin_system=create(sys,inter)`.
- Lines 29-30: Anharmonicity part from the declared interactions; implemented by `spin_system=assume(spin_system,'cavity')`.
- Lines 33-34: Transmon operators; implemented by `Cr=operator(spin_system,'C',1)`.
- Lines 38-39: Free-evolution parameters; implemented by `detunings=linspace(-20e6,20e6,256)`.
- Lines 42-43: Nominal pi/2 pulse propagator; implemented by `X=(Cr+An)/2; U90=expm(-1i*(pi/2)*full(X))`.
- Lines 45-46: Initial density matrix; implemented by `rho=state(spin_system,'BL1',1)`.
- Lines 49-50: Detection state; implemented by `L2=state(spin_system,'BL2',1)`.
- Lines 52-54: Loop over offsets; implemented by `answer=zeros(numel(detunings), numel(time_axis))`.
- Lines 57-58: Build the evolution Hamiltonian; implemented by `H=H0+2*pi*detunings(n)*N; H=(H+H')/2`.
- Lines 60-61: Propagate the Ramsey interferometer; implemented by `for k=1:numel(time_axis)`.
- Lines 63-64: Evolution propagator (small); implemented by `U=expm(-1i*full(H)*time_axis(k))`.
- Lines 66-67: Time evolution; implemented by `rho_t=U*rho*U'`.
- Lines 69-70: Final pulse; implemented by `rho_t=U90*rho_t*U90'`.
- Lines 72-73: Detection; implemented by `answer(n,k)=real(hdot(L2,rho_t))`.

### Control flow inferred from the code

- Line 55: `for` loop over `n=1:numel(detunings)`.
- Line 61: `for` loop over `k=1:numel(time_axis)`.

### Key state/data transformations

- Lines 12: computes `sys.magnet` using `sys.magnet=0`.
- Lines 15: computes `sys.isotopes` using `sys.isotopes={'T3'}`.
- Lines 18: computes `inter.modes.frqs` using `inter.modes.frqs={0}`.
- Lines 19: computes `inter.modes.anharms` using `inter.modes.anharms={-260e6}`.
- Lines 22: computes `bas.formalism` using `bas.formalism='zeeman-hilb'`.
- Lines 23: computes `bas.approximation` using `bas.approximation='none'`.
- Lines 26: computes `spin_system` using `spin_system=create(sys,inter)`.
- Lines 31: computes `H0` using `H0=hamiltonian(spin_system)`.
- Lines 34: computes `Cr` using `Cr=operator(spin_system,'C',1)`.
- Lines 35: computes `An` using `An=operator(spin_system,'A',1)`.
- Lines 36: computes `N` using `N=operator(spin_system,'N',1)`.
- Lines 39: computes `detunings` using `detunings=linspace(-20e6,20e6,256)`.
- Lines 40: computes `time_axis` using `time_axis=linspace(0,1.0e-6,256)`.
- Lines 43: computes `X` using `X=(Cr+An)/2; U90=expm(-1i*(pi/2)*full(X))`.
- Lines 46: computes `rho` using `rho=state(spin_system,'BL1',1)`.
- Lines 50: computes `L2` using `L2=state(spin_system,'BL2',1)`.
- Lines 53-54: computes `answer` using `answer=zeros(numel(detunings), numel(time_axis))`.
- Lines 58: computes `H` using `H=H0+2*pi*detunings(n)*N; H=(H+H')/2`.

## Implementation structure

- Ramsey chevron of a three-level transmon in the Duffing ap-
- proximation. A nominal pi/2 pulse prepares a coherence, and
- detuning during free evolution produces Ramsey fringes.
- Calculation time: seconds
- Magnet field
- Particle specification
- Transmon in the rotating frame
- Formalism and basis
- Spinach housekeeping
- Anharmonicity part from the declared interactions
- Transmon operators
- Free-evolution parameters

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `create()`, `basis()`, `assume()`, `hamiltonian()`, `operator()`, `state()`, `detunings()`, `time_axis()`, `answer()`, `hdot()`, `kfigure()`, `kxlabel()`, `kylabel()`, `ktitle()`.
