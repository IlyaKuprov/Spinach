# examples/quantum_tech/jaynes_cummings_a.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/examples/quantum_tech/jaynes_cummings_a.m`
- Signature: `jaynes_cummings_a()`
- Total lines: 70

## Purpose

Jaynes-Cummings coupling between a spin and an electromagnetic cavity mode with five population numbers included. The avoided crossing in the one-photon energy level splitting of the mode is plotted with respect to the detuning. Calculation time: seconds

## Physical / mathematical content

- Quantum-technology examples. The files in this area model cavity QED, transmon qubits, NV centres, and related effective Hamiltonians. The recurring mathematics is finite-dimensional quantum dynamics with ladder operators, rotating-wave-style couplings, anharmonic oscillator terms, avoided crossings, and coherent control in coupled few-mode systems.
- The physics is Jaynes-Cummings-like cavity QED: a two-level or few-level matter degree of freedom exchanges excitations with a quantised harmonic mode through rotating terms such as a†σ_- + aσ_+.

## Numerical / algorithmic content

- An eigenvalue problem is solved or analysed, so the file is extracting spectra, stationary states, avoided crossings, or modal structure from the effective Hamiltonian or superoperator.

## Code-derived implementation details

### Comment-guided execution stages

- Lines 12-13: Magnet field; implemented by `sys.magnet=0.33`.
- Lines 15-16: System; implemented by `sys.isotopes={'E','C5'}`.
- Lines 18-19: Cavity resonant with the electron; implemented by `e_frq=-sys.magnet*spin('E')/(2*pi)`.
- Lines 24-25: Basis set; implemented by `bas.formalism='zeeman-hilb'`.
- Lines 28-29: Spinach housekeeping; implemented by `spin_system=create(sys,inter)`.
- Lines 32-33: Rotating frame Hamiltonian; implemented by `spin_system=assume(spin_system,'cavity')`.
- Lines 36-37: Electron detuning operator; implemented by `Ez=operator(spin_system,{'Lz'},{1})`.
- Lines 39-40: Locate the one-excitation manifold; implemented by `spin_exc=state(spin_system,{'ZL2','BL1'},{1,2})`.
- Lines 46-47: Detuning range; implemented by `delta=2*pi*linspace(-15e6,15e6,100)`.
- Lines 49-50: Eigenvalue array; implemented by `eig_array=zeros(2,100)`.
- Lines 52-53: Loop over detunings; implemented by `for n=1:numel(delta)`.
- Lines 55-56: Make the Hamiltonian; implemented by `H=delta(n)*Ez+H_JC; H=(H+H')/2`.
- Lines 58-59: Diagonalise the one-photon manifold; implemented by `eig_array(:,n)=eig(full(one_quant'*H*one_quant))`.
- Lines 63-65: Plot one-photon case; implemented by `kfigure(); plot(1e-6*delta/(2*pi), 1e-6*eig_array/(2*pi))`.

### Control flow inferred from the code

- Line 53: `for` loop over `n=1:numel(delta)`.

### Key state/data transformations

- Lines 13: computes `sys.magnet` using `sys.magnet=0.33`.
- Lines 16: computes `sys.isotopes` using `sys.isotopes={'E','C5'}`.
- Lines 19: computes `e_frq` using `e_frq=-sys.magnet*spin('E')/(2*pi)`.
- Lines 20: computes `inter.modes.frqs` using `inter.modes.frqs={[] e_frq}`.
- Lines 21: computes `inter.modes.exchange` using `inter.modes.exchange=cell(2,2)`.
- Lines 22: computes `inter.modes.exchange{1,2}` using `inter.modes.exchange{1,2}=2.828e6`.
- Lines 25: computes `bas.formalism` using `bas.formalism='zeeman-hilb'`.
- Lines 26: computes `bas.approximation` using `bas.approximation='none'`.
- Lines 29: computes `spin_system` using `spin_system=create(sys,inter)`.
- Lines 34: computes `H_JC` using `H_JC=hamiltonian(spin_system)`.
- Lines 37: computes `Ez` using `Ez=operator(spin_system,{'Lz'},{1})`.
- Lines 40: computes `spin_exc` using `spin_exc=state(spin_system,{'ZL2','BL1'},{1,2})`.
- Lines 41: computes `cav_exc` using `cav_exc=state(spin_system,{'ZL1','BL2'},{1,2})`.
- Lines 42: computes `one_quant` using `one_quant=speye(size(H_JC,1))`.
- Lines 47: computes `delta` using `delta=2*pi*linspace(-15e6,15e6,100)`.
- Lines 50: computes `eig_array` using `eig_array=zeros(2,100)`.
- Lines 56: computes `H` using `H=delta(n)*Ez+H_JC; H=(H+H')/2`.
- Lines 59: computes `eig_array(:,n)` using `eig_array(:,n)=eig(full(one_quant'*H*one_quant))`.

## Implementation structure

- Jaynes-Cummings coupling between a spin and an electromagnetic
- cavity mode with five population numbers included. The avoided
- crossing in the one-photon energy level splitting of the mode
- is plotted with respect to the detuning.
- Calculation time: seconds
- Magnet field
- System
- Cavity resonant with the electron
- Basis set
- Spinach housekeeping
- Rotating frame Hamiltonian
- Electron detuning operator

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `spin()`, `create()`, `basis()`, `assume()`, `hamiltonian()`, `operator()`, `state()`, `speye()`, `one_quant()`, `delta()`, `eig_array()`, `kfigure()`, `kxlabel()`, `kylabel()`.
