# examples/giant_spin/case_studies/ho_pzdo4_powder.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/examples/giant_spin/case_studies/ho_pzdo4_powder.m`
- Signature: `ho_pzdo4_powder()`
- Total lines: 102

## Purpose

Powder-averaged pulsed-field magnetisation of the Ho(pzdo)4 metal- organic framework, a J=8 giant spin with a crystal field to twelfth spherical rank, under a 10 T/ms linear sweep at 2 K with spin-phonon relaxation in the generalised Lindblad form of Saito and Miyashita. The magnetisation of every orientation of the two-angle Lebedev grid is plotted alongside the powder average and the thermal equilibrium magnetisati

## Physical / mathematical content

- Giant-spin examples. The effective model treats lanthanides or high-spin centres using crystal-field / Stevens-operator Hamiltonians, Zeeman splitting, and magnetisation dynamics.

## Numerical / algorithmic content

- An eigenvalue problem is solved or analysed, so the file is extracting spectra, stationary states, avoided crossings, or modal structure from the effective Hamiltonian or superoperator.

## Code-derived implementation details

### Comment-guided execution stages

- Lines 17-18: Crystal field parameters, cm^-1, ranks 2 to 12 in Stevens operator convention; implemented by `[ks,qs,bkq]=ho_pzdo4_params()`.
- Lines 20-21: Convert Stevens coefficients into spherical tensor coefficients, Hz, rank by rank; implemented by `coeff=cell(1,12); euler=cell(1,12)`.
- Lines 27-28: Magnet must be 1 Tesla, the field is set by the sweep; implemented by `sys.magnet=1.0`.
- Lines 30-31: Parallel pool size; implemented by `sys.parallel={'processes',64}`.
- Lines 33-34: J=8 giant spin, effective g-factor 1.24, phonon bath at 2 K; implemented by `sys.isotopes={'E17'}`.
- Lines 40-41: Formalism and basis set; implemented by `bas.formalism='zeeman-hilb'`.
- Lines 44-45: Spinach housekeeping; implemented by `spin_system=create(sys,inter)`.
- Lines 48-49: Spin-phonon coupling operator: unit elements between adjacent m_J states; implemented by `Jz=full(stevens(17,1,0)); mj=diag(Jz)`.
- Lines 52-53: Super-Ohmic bath, lambda^2*I0 of the paper (lambda=10 cm^-1, I0=1e-14 ps/rad) in rad/s units; implemented by `parameters.phonon_alpha=2`.
- Lines 56-57: Observable: magnetic moment along Z in Bohr magnetons; implemented by `parameters.coil=-1.24*Jz`.
- Lines 59-60: Sweep: 10 T/ms to 10 T in 10 ns stairs, output every 1000 stairs; implemented by `parameters.field_prof=@(t) 1e4*t`.
- Lines 63-64: Two-angle Lebedev grid, outputs of every orientation returned separately; implemented by `parameters.spins={'E17'}; parameters.grid='leb_2ang_rank_29'`.
- Lines 67-68: Run the simulation; implemented by `[answers,sph_grid]=powder(spin_system,@pulsed_field,parameters,'labframe')`.
- Lines 70-71: Powder average of the magnetisation; implemented by `fields=answers{1}.field; m_avg=zeros(size(fields))`.
- Lines 76-77: Thermal equilibrium magnetisation, powder averaged, at the same fields; implemented by `[I,Q]=hamiltonian(assume(spin_system,'labframe')); [ZI,ZQ]=hamiltonian(assume(spin_system,'labframe','zeeman'))`.
- Lines 89-90: Plot the single orientation curves, the powder average, and the equilibrium; implemented by `kfigure(); hold on`.
- Lines 98-99: Save the curves; implemented by `save('ho_pzdo4_powder.mat','fields','m_avg','m_eq','answers','sph_grid')`.

### Control flow inferred from the code

- Line 22: `for` loop over `k=1:12`.
- Line 72: `for` loop over `n=1:numel(answers)`.
- Line 79: `for` loop over `n=1:numel(answers)`.
- Line 82: `for` loop over `k=1:numel(fields)`.
- Line 91: `for` loop over `n=1:numel(answers)`.

### Key state/data transformations

- Lines 18: computes `[ks,qs,bkq]` using `[ks,qs,bkq]=ho_pzdo4_params()`.
- Lines 21: computes `coeff` using `coeff=cell(1,12); euler=cell(1,12)`.
- Lines 24: computes `coeff{k}` using `coeff{k}=stev2sph(k,icm2hz(stev)); euler{k}=[0 0 0]`.
- Lines 28: computes `sys.magnet` using `sys.magnet=1.0`.
- Lines 31: computes `sys.parallel` using `sys.parallel={'processes',64}`.
- Lines 34: computes `sys.isotopes` using `sys.isotopes={'E17'}`.
- Lines 35: computes `inter.zeeman.scalar` using `inter.zeeman.scalar={1.24}`.
- Lines 36: computes `inter.giant.coeff` using `inter.giant.coeff={coeff}`.
- Lines 37: computes `inter.giant.euler` using `inter.giant.euler={euler}`.
- Lines 38: computes `inter.temperature` using `inter.temperature=2.0`.
- Lines 41: computes `bas.formalism` using `bas.formalism='zeeman-hilb'`.
- Lines 42: computes `bas.approximation` using `bas.approximation='none'`.
- Lines 45: computes `spin_system` using `spin_system=create(sys,inter)`.
- Lines 49: computes `Jz` using `Jz=full(stevens(17,1,0)); mj=diag(Jz)`.
- Lines 53: computes `parameters.phonon_alpha` using `parameters.phonon_alpha=2`.
- Lines 54: computes `parameters.phonon_i0` using `parameters.phonon_i0=1e2*1e-14*1e12*(1e-12)^2*0.1883651568463003^2`.
- Lines 57: computes `parameters.coil` using `parameters.coil=-1.24*Jz`.
- Lines 60: computes `parameters.field_prof` using `parameters.field_prof=@(t) 1e4*t`.

## Implementation structure

- Powder-averaged pulsed-field magnetisation of the Ho(pzdo)4 metal-
- organic framework, a J=8 giant spin with a crystal field to twelfth
- spherical rank, under a 10 T/ms linear sweep at 2 K with spin-phonon
- relaxation in the generalised Lindblad form of Saito and Miyashita.
- The magnetisation of every orientation of the two-angle Lebedev grid
- is plotted alongside the powder average and the thermal equilibrium
- magnetisation. Reproduces Figure 3 of
- Calculation time: hours
- Crystal field parameters, cm^-1, ranks 2 to 12 in Stevens operator convention
- Convert Stevens coefficients into spherical tensor coefficients, Hz, rank by rank
- Magnet must be 1 Tesla, the field is set by the sweep
- Parallel pool size

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `ho_pzdo4_params()`, `stev()`, `bkq()`, `stev2sph()`, `icm2hz()`, `create()`, `basis()`, `stevens()`, `double()`, `powder()`, `hamiltonian()`, `assume()`, `orientation()`, `fields()`, `m_eq()`, `kfigure()`.
