# examples/relaxation_theory/trosy_methyl.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/examples/relaxation_theory/trosy_methyl.m`
- Signature: `trosy_methyl()`
- Total lines: 155

## Purpose

Methyl trosy in a rapidly rotating 13CH3 group of a slowly tumbling protein, simulated using the Fokker-Planck formalism.

## Physical / mathematical content

- Relaxation-theory examples. The mathematical backbone is Bloch-Redfield-Wangsness or stochastic Liouville theory, spectral densities, cross-correlation terms, motional models, and extraction of longitudinal/transverse decay behaviour from superoperators.
- The file uses a Fokker-Planck-style enlarged state space in which spatial or orientational coordinates are promoted to extra dimensions and coupled to spin dynamics through differential operators.

## Numerical / algorithmic content

- The file is built around the standard Spinach workflow: create the spin system, choose a basis or context, assemble operators/superoperators, then propagate or analyse the resulting dynamics.

## Code-derived implementation details

### Comment-guided execution stages

- Lines 9-10: Magnet field; implemented by `sys.magnet=14.1`.
- Lines 12-13: Cartesian coordinates; implemented by `me_xyz=[-1.290651 1.279824 -0.248354`.
- Lines 18-19: Absolute shielding tensors (DFT); implemented by `shielding=cell(1,4)`.
- Lines 33-36: Convert shielding tensors into chemical shift tensors and put them on resonance; implemented by `shift_tensor=cell(1,4)`.
- Lines 41-42: Methyl proton chemical shifts (guess); implemented by `shift_tensor{2}=shift_tensor{2}+0.8*eye(3)`.
- Lines 46-47: J-couplings; implemented by `j_coupling=cell(4,4)`.
- Lines 55-56: Spin system instances for the three methyl rotamers; implemented by `sys.isotopes=repmat({'13C','1H','1H','1H'},1,3)`.
- Lines 61-62: Assign spins to rotamers; implemented by `inter.chem.parts={[1 2 3 4], [5 6 7 8], [9 10 11 12]}`.
- Lines 64-65: Equal rotamer populations; implemented by `inter.chem.concs=[1 1 1]`.
- Lines 67-68: Interactions within methyl rotamer A; implemented by `inter.coordinates{1}=me_xyz(1,:)`.
- Lines 78-79: Interactions within methyl rotamer B; implemented by `inter.coordinates{5}=me_xyz(1,:)`.
- Lines 89-90: Interactions within methyl rotamer C; implemented by `inter.coordinates{9}=me_xyz(1,:)`.
- Lines 100-101: Formalism and basis set; implemented by `bas.formalism='sphten-liouv'`.
- Lines 104-105: High accuracy; implemented by `sys.disable={'zte'}`.
- Lines 107-108: Methyl turning generator; implemented by `tau_m=1e-11; k_jump=1/(2*tau_m)`.
- Lines 113-114: Spinach housekeeping; implemented by `spin_system=create(sys,inter)`.
- Lines 117-118: Get a figure going; implemented by `kfigure(); scale_figure([1.5 0.8])`.
- Lines 131-132: Call frequency domain detection; implemented by `spectrum=gridfree(spin_system,@slowpass,parameters,'nmr')`.

### Control flow inferred from the code

- Line 37: `for` loop over `n=1:4`.

### Key state/data transformations

- Lines 10: computes `sys.magnet` using `sys.magnet=14.1`.
- Lines 13: computes `me_xyz` using `me_xyz=[-1.290651 1.279824 -0.248354`.
- Lines 19: computes `shielding` using `shielding=cell(1,4)`.
- Lines 20: computes `shielding{1}` using `shielding{1}=[166.0 -3.0 8.5`.
- Lines 23: computes `shielding{2}` using `shielding{2}=[ 27.8 -0.9 1.3`.
- Lines 26: computes `shielding{3}` using `shielding{3}=[ 34.6 -2.6 -1.1`.
- Lines 29: computes `shielding{4}` using `shielding{4}=[ 28.8 1.3 1.0`.
- Lines 36: computes `shift_tensor` using `shift_tensor=cell(1,4)`.
- Lines 38: computes `shift_tensor{n}` using `shift_tensor{n}=-remtrace(shielding{n})`.
- Lines 42: computes `shift_tensor{2}` using `shift_tensor{2}=shift_tensor{2}+0.8*eye(3)`.
- Lines 43: computes `shift_tensor{3}` using `shift_tensor{3}=shift_tensor{3}+1.0*eye(3)`.
- Lines 44: computes `shift_tensor{4}` using `shift_tensor{4}=shift_tensor{4}+1.2*eye(3)`.
- Lines 47: computes `j_coupling` using `j_coupling=cell(4,4)`.
- Lines 48: computes `j_coupling{1,2}` using `j_coupling{1,2}=125`.
- Lines 49: computes `j_coupling{1,3}` using `j_coupling{1,3}=125`.
- Lines 50: computes `j_coupling{1,4}` using `j_coupling{1,4}=125`.
- Lines 51: computes `j_coupling{2,3}` using `j_coupling{2,3}=-12`.
- Lines 52: computes `j_coupling{2,4}` using `j_coupling{2,4}=-12`.

## Implementation structure

- Methyl trosy in a rapidly rotating 13CH3 group
- of a slowly tumbling protein, simulated using
- the Fokker-Planck formalism.
- Magnet field
- Cartesian coordinates
- Absolute shielding tensors (DFT)
- Convert shielding tensors into
- chemical shift tensors and put
- them on resonance
- Methyl proton chemical shifts (guess)
- J-couplings
- Spin system instances for the three methyl rotamers

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `remtrace()`, `me_xyz()`, `create()`, `basis()`, `kfigure()`, `scale_figure()`, `state()`, `gridfree()`, `subplot()`, `plot_1d()`.
