# examples/giant_spin/case_studies/mn3_trimer_levels.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/examples/giant_spin/case_studies/mn3_trimer_levels.m`
- Signature: `mn3_trimer_levels()`
- Total lines: 63

## Purpose

Zeeman energy level diagram of the (CH6N3)2MnCl4 molecular crystal, a linear trimer of three S=5/2 manganese ions with isotropic exchange between neighbours and an axial plus rhombic zero-field splitting on every ion, from zero to 10 Tesla in the full 216-state Hilbert space. The lowest levels are the ones that the 16-state and 26-state effec- tive bases of the paper are built to reproduce. Reproduces Figure 5 of Cal

## Physical / mathematical content

- Giant-spin examples. The effective model treats lanthanides or high-spin centres using crystal-field / Stevens-operator Hamiltonians, Zeeman splitting, and magnetisation dynamics.

## Numerical / algorithmic content

- An eigenvalue problem is solved or analysed, so the file is extracting spectra, stationary states, avoided crossings, or modal structure from the effective Hamiltonian or superoperator.

## Code-derived implementation details

### Comment-guided execution stages

- Lines 17-18: Magnet must be 1 Tesla, the field is set below; implemented by `sys.magnet=1.0`.
- Lines 20-21: Parallel pool size; implemented by `sys.parallel={'processes',4}`.
- Lines 23-24: Three S=5/2 spins with g=2; implemented by `sys.isotopes={'E6','E6','E6'}`.
- Lines 27-28: Isotropic exchange, J=-2.42 cm^-1 in the H=-2*J*S1*S2 convention of the paper; implemented by `inter.coupling.matrix=cell(3,3)`.
- Lines 32-33: Zero-field splitting, D=0.167 cm^-1 and E=0.040 cm^-1 on every ion; implemented by `for n=1:3`.
- Lines 37-38: Formalism and basis set; implemented by `bas.formalism='zeeman-hilb'`.
- Lines 41-42: Spinach housekeeping; implemented by `spin_system=create(sys,inter)`.
- Lines 45-46: Field-free Hamiltonian and the Zeeman operator per Tesla; implemented by `[I,Q]=hamiltonian(assume(spin_system,'labframe')); H0=I+orientation(Q,[0 0 0])`.
- Lines 49-50: Energy levels on a field grid, cm^-1; implemented by `fields=linspace(0,10,201); levels=zeros(size(H0,1),numel(fields))`.
- Lines 55-56: Plot the lowest thirty levels relative to the field-free ground state; implemented by `kfigure(); plot(fields,levels(1:30,:)-levels(1,1)); kgrid; xlim tight; ylim padded`.
- Lines 59-60: Save the levels; implemented by `save('mn3_trimer_levels.mat','fields','levels')`.

### Control flow inferred from the code

- Line 33: `for` loop over `n=1:3`.
- Line 51: `for` loop over `k=1:numel(fields)`.

### Key state/data transformations

- Lines 18: computes `sys.magnet` using `sys.magnet=1.0`.
- Lines 21: computes `sys.parallel` using `sys.parallel={'processes',4}`.
- Lines 24: computes `sys.isotopes` using `sys.isotopes={'E6','E6','E6'}`.
- Lines 25: computes `inter.zeeman.scalar` using `inter.zeeman.scalar={2.0 2.0 2.0}`.
- Lines 28: computes `inter.coupling.matrix` using `inter.coupling.matrix=cell(3,3)`.
- Lines 29: computes `inter.coupling.matrix{1,2}` using `inter.coupling.matrix{1,2}=-2*icm2hz(-2.42)*eye(3)`.
- Lines 30: computes `inter.coupling.matrix{2,3}` using `inter.coupling.matrix{2,3}=-2*icm2hz(-2.42)*eye(3)`.
- Lines 34: computes `inter.coupling.matrix{n,n}` using `inter.coupling.matrix{n,n}=zfs2mat(icm2hz(0.167),icm2hz(0.040),0,0,0)`.
- Lines 38: computes `bas.formalism` using `bas.formalism='zeeman-hilb'`.
- Lines 39: computes `bas.approximation` using `bas.approximation='none'`.
- Lines 42: computes `spin_system` using `spin_system=create(sys,inter)`.
- Lines 46: computes `[I,Q]` using `[I,Q]=hamiltonian(assume(spin_system,'labframe')); H0=I+orientation(Q,[0 0 0])`.
- Lines 47: computes `Z` using `Z=hamiltonian(assume(spin_system,'labframe','zeeman')); H0=H0-Z`.
- Lines 50: computes `fields` using `fields=linspace(0,10,201); levels=zeros(size(H0,1),numel(fields))`.
- Lines 52: computes `H` using `H=full(H0+fields(k)*Z); levels(:,k)=hz2icm(sort(eig((H+H')/2))/(2*pi))`.

## Implementation structure

- Zeeman energy level diagram of the (CH6N3)2MnCl4 molecular crystal, a
- linear trimer of three S=5/2 manganese ions with isotropic exchange
- between neighbours and an axial plus rhombic zero-field splitting on
- every ion, from zero to 10 Tesla in the full 216-state Hilbert space.
- The lowest levels are the ones that the 16-state and 26-state effec-
- tive bases of the paper are built to reproduce. Reproduces Figure 5
- Calculation time: seconds
- Magnet must be 1 Tesla, the field is set below
- Parallel pool size
- Three S=5/2 spins with g=2
- Isotropic exchange, J=-2.42 cm^-1 in the H=-2*J*S1*S2 convention of the paper
- Zero-field splitting, D=0.167 cm^-1 and E=0.040 cm^-1 on every ion

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `icm2hz()`, `zfs2mat()`, `create()`, `basis()`, `hamiltonian()`, `assume()`, `orientation()`, `fields()`, `levels()`, `hz2icm()`, `kfigure()`, `kxlabel()`, `kylabel()`, `save()`.
