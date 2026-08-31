# examples/nmr_solids/case_studies/cp_square_vs_ramp/cp_square_vs_ramp.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/examples/nmr_solids/case_studies/cp_square_vs_ramp/cp_square_vs_ramp.m`
- Signature: `cp_square_vs_ramp()`
- Total lines: 81

## Purpose

1H-15N cross-polarisation experiment in the doubly rotating frame using (a) fixed amplitude CP; (b) linearly ramped CP; (c) tangent-ramped CP. Static powder simulation demonstra- ting the advantages of ramped cross-polarisation. For fur- ther information, see: Calculation time: seconds

## Physical / mathematical content

- Solid-state NMR examples. The key physics is anisotropic spin interactions under static or magic-angle-spinning conditions: chemical-shift anisotropy, dipolar coupling, quadrupolar coupling, cross-polarisation, and orientation averaging using Floquet, Fokker-Planck, or direct powder quadrature formalisms.

## Numerical / algorithmic content

- The file is built around the standard Spinach workflow: create the spin system, choose a basis or context, assemble operators/superoperators, then propagate or analyse the resulting dynamics.

## Code-derived implementation details

### Comment-guided execution stages

- Lines 15-16: System specification; implemented by `sys.magnet=9.394`.
- Lines 19-20: Interactions; implemented by `inter.zeeman.scalar={0.00 0.00}`.
- Lines 25-26: Basis set; implemented by `bas.formalism='sphten-liouv'`.
- Lines 29-30: Spinach housekeeping; implemented by `spin_system=create(sys,inter)`.
- Lines 33-34: Common experiment parameters; implemented by `parameters.time_steps=2e-6*ones(1,500)`.
- Lines 43-44: Simulate fixed amplitude CP; implemented by `irr_powers_a=[5e4*ones(1,500); 5e4*ones(1,500)]`.
- Lines 48-49: Simulate linearly ramped amplitude CP; implemented by `ramp_up=linspace(0,1,500); ramp_down=fliplr(ramp_up)`.
- Lines 54-55: Simulate tangent ramped amplitude CP; implemented by `ramp_up=tan(linspace(-1.4,1.4,500))`.
- Lines 63-64: Plotting; implemented by `kfigure(); scale_figure([1.5 0.75])`.

### Key state/data transformations

- Lines 16: computes `sys.magnet` using `sys.magnet=9.394`.
- Lines 17: computes `sys.isotopes` using `sys.isotopes={'15N','1H'}`.
- Lines 20: computes `inter.zeeman.scalar` using `inter.zeeman.scalar={0.00 0.00}`.
- Lines 21: computes `inter.coordinates` using `inter.coordinates={[0.00 0.00 0.00]`.
- Lines 23: computes `inter.temperature` using `inter.temperature=298`.
- Lines 26: computes `bas.formalism` using `bas.formalism='sphten-liouv'`.
- Lines 27: computes `bas.approximation` using `bas.approximation='none'`.
- Lines 30: computes `spin_system` using `spin_system=create(sys,inter)`.
- Lines 34: computes `parameters.time_steps` using `parameters.time_steps=2e-6*ones(1,500)`.
- Lines 35-36: computes `parameters.irr_opers` using `parameters.irr_opers={operator(spin_system,'Ly','1H') operator(spin_system,'Lx','15N')}`.
- Lines 37: computes `parameters.exc_opers` using `parameters.exc_opers={operator(spin_system,'Lx','1H')}`.
- Lines 38: computes `parameters.coil` using `parameters.coil=state(spin_system,'Lx','15N')`.
- Lines 39: computes `parameters.grid` using `parameters.grid='rep_2ang_6400pts_sph'`.
- Lines 40: computes `parameters.needs` using `parameters.needs={'aniso_eq'}`.
- Lines 41: computes `parameters.spins` using `parameters.spins={'15N'}`.
- Lines 44: computes `irr_powers_a` using `irr_powers_a=[5e4*ones(1,500); 5e4*ones(1,500)]`.
- Lines 45: computes `parameters.irr_powers` using `parameters.irr_powers=irr_powers_a`.
- Lines 46: computes `fid_a` using `fid_a=powder(spin_system,@cp_contact_hard,parameters,'nmr')`.

## Implementation structure

- 1H-15N cross-polarisation experiment in the doubly rotating
- frame using (a) fixed amplitude CP; (b) linearly ramped CP;
- (c) tangent-ramped CP. Static powder simulation demonstra-
- ting the advantages of ramped cross-polarisation. For fur-
- ther information, see:
- Calculation time: seconds
- System specification
- Interactions
- Basis set
- Spinach housekeeping
- Common experiment parameters
- Simulate fixed amplitude CP

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `create()`, `basis()`, `operator()`, `state()`, `powder()`, `fliplr()`, `kfigure()`, `scale_figure()`, `cumsum()`, `subplot()`, `time_axis()`, `kxlabel()`, `kylabel()`, `ylim()`, `klegend()`.
