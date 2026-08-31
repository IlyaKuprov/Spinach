# examples/quantum_tech/tavis_cummings_splitting.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/examples/quantum_tech/tavis_cummings_splitting.m`
- Signature: `tavis_cummings_splitting()`
- Total lines: 86

## Purpose

Collective normal-mode splitting in the Tavis-Cummings model for one to four identical electron spins coupled to a common microwave cavity mode. The bright-state splitting follows the square-root scaling of Tavis and Cummings, Phys. Rev. 170, 379 (1968). Calculation time: seconds

## Physical / mathematical content

- Quantum-technology examples. The files in this area model cavity QED, transmon qubits, NV centres, and related effective Hamiltonians. The recurring mathematics is finite-dimensional quantum dynamics with ladder operators, rotating-wave-style couplings, anharmonic oscillator terms, avoided crossings, and coherent control in coupled few-mode systems.

## Numerical / algorithmic content

- An eigenvalue problem is solved or analysed, so the file is extracting spectra, stationary states, avoided crossings, or modal structure from the effective Hamiltonian or superoperator.

## Code-derived implementation details

### Comment-guided execution stages

- Lines 13-14: Coupling strength; implemented by `coupling=3e6`.
- Lines 16-17: Preallocate the splitting array; implemented by `spin_counts=1:4`.
- Lines 20-21: Loop over ensemble sizes; implemented by `for n=spin_counts`.
- Lines 23-24: Magnet field; implemented by `sys.magnet=0`.
- Lines 26-27: Particle specification; implemented by `sys.isotopes=[repmat({'E'},1,n) {'C3'}]`.
- Lines 29-30: Resonant cavity coupled to every spin; implemented by `clear('inter')`.
- Lines 37-38: Formalism and basis; implemented by `bas.formalism='zeeman-hilb'`.
- Lines 41-42: Spinach housekeeping; implemented by `spin_system=create(sys,inter)`.
- Lines 45-46: Tavis-Cummings Hamiltonian from the declared couplings; implemented by `spin_system=assume(spin_system,'cavity')`.
- Lines 50-51: Locate the one-excitation manifold; implemented by `one_quant=speye(size(H,1))`.
- Lines 63-64: Extract the bright-mode splitting; implemented by `energies=sort(real(eig(full(one_quant'*H*one_quant))))`.
- Lines 69-70: Compute the analytical splitting; implemented by `analytical=2*sqrt(spin_counts)*2*pi*coupling`.
- Lines 72-73: Validate the square-root scaling; implemented by `if max(abs(splitting-analytical)./analytical)>1e-10`.
- Lines 77-78: Plot numerical and analytical splittings; implemented by `kfigure(); plot(spin_counts,splitting/(2*pi*1e6),'ko','MarkerFaceColor','k')`.

### Control flow inferred from the code

- Line 21: `for` loop over `n=spin_counts`.
- Line 33: `for` loop over `k=1:n`.
- Line 53: `for` loop over `k=1:n`.
- Line 73: conditional branch on `max(abs(splitting-analytical)./analytical)>1e-10`.

### Key state/data transformations

- Lines 14: computes `coupling` using `coupling=3e6`.
- Lines 17: computes `spin_counts` using `spin_counts=1:4`.
- Lines 18: computes `splitting` using `splitting=zeros(size(spin_counts))`.
- Lines 24: computes `sys.magnet` using `sys.magnet=0`.
- Lines 27: computes `sys.isotopes` using `sys.isotopes=[repmat({'E'},1,n) {'C3'}]`.
- Lines 31: computes `inter.modes.frqs` using `inter.modes.frqs=[repmat({[]},1,n) {0}]`.
- Lines 32: computes `inter.modes.exchange` using `inter.modes.exchange=cell(n+1,n+1)`.
- Lines 34: computes `inter.modes.exchange{k,n+1}` using `inter.modes.exchange{k,n+1}=coupling`.
- Lines 38: computes `bas.formalism` using `bas.formalism='zeeman-hilb'`.
- Lines 39: computes `bas.approximation` using `bas.approximation='none'`.
- Lines 42: computes `spin_system` using `spin_system=create(sys,inter)`.
- Lines 47: computes `H` using `H=hamiltonian(spin_system)`.
- Lines 51: computes `one_quant` using `one_quant=speye(size(H,1))`.
- Lines 52: computes `subspace` using `subspace=zeros(n+1,1)`.
- Lines 54: computes `spin_state` using `spin_state=repmat({'ZL1'},1,n)`.
- Lines 55: computes `spin_state{k}` using `spin_state{k}='ZL2'`.
- Lines 56: computes `projector` using `projector=state(spin_system,[spin_state {'BL1'}],num2cell(1:(n+1)))`.
- Lines 57: computes `subspace(k)` using `subspace(k)=find(diag(projector)>0.5)`.

## Implementation structure

- Collective normal-mode splitting in the Tavis-Cummings model
- for one to four identical electron spins coupled to a common
- microwave cavity mode. The bright-state splitting follows the
- square-root scaling of Tavis and Cummings, Phys. Rev. 170,
- 379 (1968).
- Calculation time: seconds
- Coupling strength
- Preallocate the splitting array
- Loop over ensemble sizes
- Magnet field
- Particle specification
- Resonant cavity coupled to every spin

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `clear()`, `create()`, `basis()`, `assume()`, `hamiltonian()`, `speye()`, `state()`, `num2cell()`, `subspace()`, `one_quant()`, `splitting()`, `kfigure()`, `kxlabel()`, `kylabel()`, `ktitle()`, `klegend()`.
