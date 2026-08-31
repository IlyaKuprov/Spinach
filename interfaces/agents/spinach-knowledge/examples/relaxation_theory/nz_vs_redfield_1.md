# examples/relaxation_theory/nz_vs_redfield_1.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/examples/relaxation_theory/nz_vs_redfield_1.m`
- Signature: `nz_vs_redfield_1()`
- Total lines: 98

## Purpose

Nakajima-Zwanzig relaxation theory against Redfield theory for a two-spin system with dipolar and CSA cross-correlations. The three superoperators compared are the off-shell NZ kernel (resolvent form), the on-shell NZ kernel (back-rotated form), and Redfield theory. The on-shell kernel at zero shift reproduces Redfield theory exactly; the off-shell kernel agrees with Redfield theory on the zero-frequency subspace of 

## Physical / mathematical content

- Relaxation-theory examples. The mathematical backbone is Bloch-Redfield-Wangsness or stochastic Liouville theory, spectral densities, cross-correlation terms, motional models, and extraction of longitudinal/transverse decay behaviour from superoperators.
- The relaxation model is Redfield-type perturbation theory: fluctuating interactions enter through correlation functions or spectral densities and generate a linear relaxation superoperator.
- The spin physics includes through-space magnetic dipole-dipole coupling, a rank-2 anisotropic interaction with strong orientation dependence and characteristic secular/non-secular structure.
- Chemical-shift anisotropy is present: shielding is treated as a second-rank tensor whose orientation relative to the field or rotor axis modulates line shapes and transfer dynamics.

## Numerical / algorithmic content

- The file is built around the standard Spinach workflow: create the spin system, choose a basis or context, assemble operators/superoperators, then propagate or analyse the resulting dynamics.

## Code-derived implementation details

### Comment-guided execution stages

- Lines 18-19: Magnet and isotopes; implemented by `sys.magnet=14.1`.
- Lines 22-23: Chemical shielding tensors; implemented by `inter.zeeman.eigs={[7 15 -22]`.
- Lines 28-29: Coordinates for the dipolar coupling; implemented by `inter.coordinates={[0.00 0.00 0.00]`.
- Lines 32-33: Relaxation theory common settings; implemented by `inter.equilibrium='zero'`.
- Lines 37-38: Basis set; implemented by `bas.formalism='sphten-liouv'`.
- Lines 41-42: Redfield superoperator; implemented by `inter_rf=inter; inter_rf.relaxation={'redfield'}; inter_rf.tau_c={200e-12}`.
- Lines 45-46: On-shell NZ superoperator at zero shift; implemented by `inter_on=inter; inter_on.relaxation={'naka-zwan'}; inter_on.tau_c={200e-12}`.
- Lines 50-51: Off-shell NZ superoperator at zero shift; implemented by `inter_off=inter_on; inter_off.nz_onshell=false`.
- Lines 55-57: On-shell kernel at zero shift is Redfield theory; implemented by `disp(['on-shell NZ vs Redfield, relative difference: ' num2str(norm(R_on-R_rf,1)/norm(R_rf,1))])`.
- Lines 59-60: Off-shell kernel agrees with Redfield theory on the kite; implemented by `V0=null(full(hamiltonian(assume(spin_system,'labframe'))))`.
- Lines 66-67: Off-shell difference is first order in omega*tau_c; implemented by `tau_grid=[2.5e-12 5e-12 10e-12 20e-12]; rel_diff=zeros(1,4)`.
- Lines 78-79: Lifetime shift suppresses relaxation rates; implemented by `shift_grid=[0 1e9 1e10 1e11]; max_rates=zeros(1,4)`.
- Lines 88-89: Plot the two trends; implemented by `kfigure(); scale_figure([1.75 0.75])`.

### Control flow inferred from the code

- Line 68: `for` loop over `n=1:4`.
- Line 80: `for` loop over `n=1:4`.

### Key state/data transformations

- Lines 19: computes `sys.magnet` using `sys.magnet=14.1`.
- Lines 20: computes `sys.isotopes` using `sys.isotopes={'1H','13C'}`.
- Lines 23: computes `inter.zeeman.eigs` using `inter.zeeman.eigs={[7 15 -22]`.
- Lines 25: computes `inter.zeeman.euler` using `inter.zeeman.euler={[pi/3 pi/4 pi/5]`.
- Lines 29: computes `inter.coordinates` using `inter.coordinates={[0.00 0.00 0.00]`.
- Lines 33: computes `inter.equilibrium` using `inter.equilibrium='zero'`.
- Lines 34: computes `inter.rlx_keep` using `inter.rlx_keep='labframe'`.
- Lines 35: computes `inter.rlx_dfs` using `inter.rlx_dfs='keep'`.
- Lines 38: computes `bas.formalism` using `bas.formalism='sphten-liouv'`.
- Lines 39: computes `bas.approximation` using `bas.approximation='none'`.
- Lines 42: computes `inter_rf` using `inter_rf=inter; inter_rf.relaxation={'redfield'}; inter_rf.tau_c={200e-12}`.
- Lines 43: computes `R_rf` using `R_rf=full(relaxation(basis(create(sys,inter_rf),bas)))`.
- Lines 46: computes `inter_on` using `inter_on=inter; inter_on.relaxation={'naka-zwan'}; inter_on.tau_c={200e-12}`.
- Lines 47: computes `inter_on.nz_onshell` using `inter_on.nz_onshell=true; inter_on.nz_shift=0`.
- Lines 48: computes `R_on` using `R_on=full(relaxation(basis(create(sys,inter_on),bas)))`.
- Lines 51: computes `inter_off` using `inter_off=inter_on; inter_off.nz_onshell=false`.
- Lines 52: computes `spin_system` using `spin_system=basis(create(sys,inter_off),bas)`.
- Lines 53: computes `R_off` using `R_off=full(relaxation(spin_system))`.

## Implementation structure

- Nakajima-Zwanzig relaxation theory against Redfield theory for a
- two-spin system with dipolar and CSA cross-correlations. The three
- superoperators compared are the off-shell NZ kernel (resolvent form),
- the on-shell NZ kernel (back-rotated form), and Redfield theory. The
- on-shell kernel at zero shift reproduces Redfield theory exactly; the
- off-shell kernel agrees with Redfield theory on the zero-frequency
- subspace of the coherent Liouvillian and differs in first order in
- omega*tau_c on coherences; a lifetime shift suppresses all rates by
- pushing the kernel off the real axis.
- Calculation time: minutes
- Magnet and isotopes
- Chemical shielding tensors

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `relaxation()`, `basis()`, `create()`, `num2str()`, `null()`, `hamiltonian()`, `assume()`, `tau_grid()`, `rel_diff()`, `shift_grid()`, `max_rates()`, `kfigure()`, `scale_figure()`, `subplot()`, `kxlabel()`, `kylabel()`.
