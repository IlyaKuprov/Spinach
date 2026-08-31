# examples/relaxation_theory/hfc_antisymm_1.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/examples/relaxation_theory/hfc_antisymm_1.m`
- Signature: `hfc_antisymm_1()`
- Total lines: 73

## Purpose

Longitudinal and transverse relaxation rates in a system with a significant antisymmetry in the hyperfine tensor. Calculation time: seconds

## Physical / mathematical content

- Relaxation-theory examples. The mathematical backbone is Bloch-Redfield-Wangsness or stochastic Liouville theory, spectral densities, cross-correlation terms, motional models, and extraction of longitudinal/transverse decay behaviour from superoperators.
- The relaxation model is Redfield-type perturbation theory: fluctuating interactions enter through correlation functions or spectral densities and generate a linear relaxation superoperator.

## Numerical / algorithmic content

- The file is built around the standard Spinach workflow: create the spin system, choose a basis or context, assemble operators/superoperators, then propagate or analyse the resulting dynamics.

## Code-derived implementation details

### Comment-guided execution stages

- Lines 10-11: System specification; implemented by `sys.magnet=0.33`.
- Lines 18-19: Relaxation theory parameters; implemented by `inter.relaxation={'redfield'}`.
- Lines 24-25: Basis set; implemented by `bas.formalism='sphten-liouv'`.
- Lines 28-29: Spinach housekeeping; implemented by `spin_system=create(sys,inter)`.
- Lines 32-33: Relaxation superoperator; implemented by `R=relaxation(spin_system)`.
- Lines 35-37: Textbook rates; implemented by `[r1,r2,rx]=rlx_hfc(sys.magnet,2*pi*inter.coupling.matrix{1,2}, sys.isotopes,inter.tau_c{1})`.
- Lines 39-40: Textbook and Spinach R1 for first spin; implemented by `rho=state(spin_system,'Lz',1); rho=rho/norm(rho,2)`.
- Lines 45-46: Textbook and Spinach R1 for second spin; implemented by `rho=state(spin_system,'Lz',2); rho=rho/norm(rho,2)`.
- Lines 51-52: Textbook and Spinach R2 for first spin; implemented by `rho=state(spin_system,'L+',1); rho=rho/norm(rho,2)`.
- Lines 57-58: Textbook and Spinach R2 for second spin; implemented by `rho=state(spin_system,'L+',2); rho=rho/norm(rho,2)`.
- Lines 63-64: Print textbook cross-relaxation rate; implemented by `rho_a=state(spin_system,'Lz',1); rho_a=rho_a/norm(rho_a,2)`.
- Lines 69-70: Print the full relaxation superoperator; implemented by `disp('Complete relaxation superoperator, IST basis:'); disp(full(R))`.

### Key state/data transformations

- Lines 11: computes `sys.magnet` using `sys.magnet=0.33`.
- Lines 12: computes `sys.isotopes` using `sys.isotopes={'1H','E'}`.
- Lines 13: computes `inter.coupling.matrix` using `inter.coupling.matrix=cell(2,2)`.
- Lines 14: computes `inter.coupling.matrix{1,2}` using `inter.coupling.matrix{1,2}=1e6*[10.0 1.0 1.5`.
- Lines 19: computes `inter.relaxation` using `inter.relaxation={'redfield'}`.
- Lines 20: computes `inter.equilibrium` using `inter.equilibrium='zero'`.
- Lines 21: computes `inter.rlx_keep` using `inter.rlx_keep='labframe'`.
- Lines 22: computes `inter.tau_c` using `inter.tau_c={10e-12}`.
- Lines 25: computes `bas.formalism` using `bas.formalism='sphten-liouv'`.
- Lines 26: computes `bas.approximation` using `bas.approximation='none'`.
- Lines 29: computes `spin_system` using `spin_system=create(sys,inter)`.
- Lines 33: computes `R` using `R=relaxation(spin_system)`.
- Lines 36-37: computes `[r1,r2,rx]` using `[r1,r2,rx]=rlx_hfc(sys.magnet,2*pi*inter.coupling.matrix{1,2}, sys.isotopes,inter.tau_c{1})`.
- Lines 40: computes `rho` using `rho=state(spin_system,'Lz',1); rho=rho/norm(rho,2)`.
- Lines 64: computes `rho_a` using `rho_a=state(spin_system,'Lz',1); rho_a=rho_a/norm(rho_a,2)`.
- Lines 65: computes `rho_b` using `rho_b=state(spin_system,'Lz',2); rho_b=rho_b/norm(rho_b,2)`.

## Implementation structure

- Longitudinal and transverse relaxation rates in a system
- with a significant antisymmetry in the hyperfine tensor.
- Calculation time: seconds
- System specification
- Relaxation theory parameters
- Basis set
- Spinach housekeeping
- Relaxation superoperator
- Textbook rates
- Textbook and Spinach R1 for first spin
- Textbook and Spinach R1 for second spin
- Textbook and Spinach R2 for first spin

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `create()`, `basis()`, `relaxation()`, `rlx_hfc()`, `state()`, `num2str()`.
