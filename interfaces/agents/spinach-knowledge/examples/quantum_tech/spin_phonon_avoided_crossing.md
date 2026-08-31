# examples/quantum_tech/spin_phonon_avoided_crossing.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/examples/quantum_tech/spin_phonon_avoided_crossing.m`
- Signature: `spin_phonon_avoided_crossing()`
- Total lines: 78

## Purpose

Avoided crossing between an electron spin transition and a quantised phonon mode in the resonant spin-phonon exchange model. The phonon is requested with the V# particle syntax. Calculation time: seconds

## Physical / mathematical content

- Quantum-technology examples. The files in this area model cavity QED, transmon qubits, NV centres, and related effective Hamiltonians. The recurring mathematics is finite-dimensional quantum dynamics with ladder operators, rotating-wave-style couplings, anharmonic oscillator terms, avoided crossings, and coherent control in coupled few-mode systems.

## Numerical / algorithmic content

- An eigenvalue problem is solved or analysed, so the file is extracting spectra, stationary states, avoided crossings, or modal structure from the effective Hamiltonian or superoperator.

## Code-derived implementation details

### Comment-guided execution stages

- Lines 11-12: Magnet field; implemented by `sys.magnet=0`.
- Lines 14-15: Particle specification; implemented by `sys.isotopes={'E','V3'}`.
- Lines 17-18: Resonant phonon mode in the rotating frame; implemented by `inter.modes.frqs={[] 0}`.
- Lines 22-23: Formalism and basis; implemented by `bas.formalism='zeeman-hilb'`.
- Lines 26-27: Spinach housekeeping; implemented by `spin_system=create(sys,inter)`.
- Lines 30-31: Exchange Hamiltonian, 'cavity' is the set that keeps spin-mode exchange; implemented by `spin_system=assume(spin_system,'cavity')`.
- Lines 35-36: Spin operator; implemented by `electron_z=operator(spin_system,'Lz',1)`.
- Lines 38-39: Coupling parameters; implemented by `g=2*pi*inter.modes.exchange{1,2}`.
- Lines 42-43: Locate the one-excitation manifold; implemented by `spin_exc=state(spin_system,{'ZL2','BL1'},{1,2})`.
- Lines 50-51: Preallocate the eigenvalue array; implemented by `levels=zeros(2,numel(detuning))`.
- Lines 53-54: Sweep the spin-phonon detuning; implemented by `for n=1:numel(detuning)`.
- Lines 56-57: Build the rotating-frame Hamiltonian; implemented by `H=detuning(n)*electron_z+action`.
- Lines 59-60: Extract the dressed one-quantum doublet; implemented by `H=one_quant'*H*one_quant`.
- Lines 65-66: Validate the resonant avoided-crossing gap; implemented by `gap=min(levels(2,:)-levels(1,:))`.
- Lines 71-72: Plot the avoided crossing; implemented by `kfigure(); plot(detuning/(2*pi*1e6),levels,'LineWidth',1.5)`.

### Control flow inferred from the code

- Line 54: `for` loop over `n=1:numel(detuning)`.
- Line 67: conditional branch on `abs(gap-2*g/(2*pi*1e6))>1e-6`.

### Key state/data transformations

- Lines 12: computes `sys.magnet` using `sys.magnet=0`.
- Lines 15: computes `sys.isotopes` using `sys.isotopes={'E','V3'}`.
- Lines 18: computes `inter.modes.frqs` using `inter.modes.frqs={[] 0}`.
- Lines 19: computes `inter.modes.exchange` using `inter.modes.exchange=cell(2,2)`.
- Lines 20: computes `inter.modes.exchange{1,2}` using `inter.modes.exchange{1,2}=4e6`.
- Lines 23: computes `bas.formalism` using `bas.formalism='zeeman-hilb'`.
- Lines 24: computes `bas.approximation` using `bas.approximation='none'`.
- Lines 27: computes `spin_system` using `spin_system=create(sys,inter)`.
- Lines 32: computes `action` using `action=hamiltonian(spin_system)`.
- Lines 36: computes `electron_z` using `electron_z=operator(spin_system,'Lz',1)`.
- Lines 39: computes `g` using `g=2*pi*inter.modes.exchange{1,2}`.
- Lines 40: computes `detuning` using `detuning=2*pi*linspace(-20e6,20e6,121)`.
- Lines 43: computes `spin_exc` using `spin_exc=state(spin_system,{'ZL2','BL1'},{1,2})`.
- Lines 44: computes `phon_exc` using `phon_exc=state(spin_system,{'ZL1','BL2'},{1,2})`.
- Lines 45: computes `one_quant` using `one_quant=speye(size(action,1))`.
- Lines 46: computes `spin_idx` using `spin_idx=find(diag(spin_exc)>0.5)`.
- Lines 47: computes `phon_idx` using `phon_idx=find(diag(phon_exc)>0.5)`.
- Lines 51: computes `levels` using `levels=zeros(2,numel(detuning))`.

## Implementation structure

- Avoided crossing between an electron spin transition and a
- quantised phonon mode in the resonant spin-phonon exchange
- model. The phonon is requested with the V# particle syntax.
- Calculation time: seconds
- Magnet field
- Particle specification
- Resonant phonon mode in the rotating frame
- Formalism and basis
- Spinach housekeeping
- Exchange Hamiltonian, 'cavity' is the set that keeps spin-mode exchange
- Spin operator
- Coupling parameters

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `create()`, `basis()`, `assume()`, `hamiltonian()`, `operator()`, `state()`, `speye()`, `one_quant()`, `detuning()`, `levels()`, `kfigure()`, `kxlabel()`, `kylabel()`, `ktitle()`.
