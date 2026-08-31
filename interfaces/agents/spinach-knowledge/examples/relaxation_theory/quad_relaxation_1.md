# examples/relaxation_theory/quad_relaxation_1.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/examples/relaxation_theory/quad_relaxation_1.m`
- Signature: `quad_relaxation_1()`
- Total lines: 51

## Purpose

14N quadrupolar relaxation in glycine in liquid state. The numerical output of Spinach is compared to the analytical equation from the textbook. Calculation time: seconds

## Physical / mathematical content

- Relaxation-theory examples. The mathematical backbone is Bloch-Redfield-Wangsness or stochastic Liouville theory, spectral densities, cross-correlation terms, motional models, and extraction of longitudinal/transverse decay behaviour from superoperators.
- The relaxation model is Redfield-type perturbation theory: fluctuating interactions enter through correlation functions or spectral densities and generate a linear relaxation superoperator.
- Quadrupolar physics is relevant: nuclei with spin > 1/2 interact with the electric field gradient tensor, introducing second-rank anisotropy, asymmetry, and overtone or MQ phenomena.

## Numerical / algorithmic content

- The file is built around the standard Spinach workflow: create the spin system, choose a basis or context, assemble operators/superoperators, then propagate or analyse the resulting dynamics.

## Code-derived implementation details

### Comment-guided execution stages

- Lines 11-12: System specification; implemented by `sys.magnet=14.1`.
- Lines 15-16: Spin quantum number and quadrupolar tensor; implemented by `[~,s_mult]=spin(sys.isotopes{1}); s_qnum=(s_mult-1)/2`.
- Lines 19-20: Relaxation theory; implemented by `inter.relaxation={'redfield'}`.
- Lines 25-26: Basis set; implemented by `bas.formalism='sphten-liouv'`.
- Lines 29-30: Spinach housekeeping; implemented by `spin_system=create(sys,inter)`.
- Lines 33-34: Relaxation superoperator; implemented by `R=relaxation(spin_system)`.
- Lines 36-37: Textbook relaxation rate expressions; implemented by `[r1,r2]=rlx_nqi(s_qnum,14.1*spin(sys.isotopes{1}),1.18e6,0.53,1e-9)`.
- Lines 39-40: States of interest; implemented by `Lz=state(spin_system,'Lz',1)`.
- Lines 44-45: Print the answers; implemented by `disp([sys.isotopes{1} ' longitudinal relaxation rate, Spinach: ' num2str(-Lz'*R*Lz)])`.

### Key state/data transformations

- Lines 12: computes `sys.magnet` using `sys.magnet=14.1`.
- Lines 13: computes `sys.isotopes` using `sys.isotopes={'14N'}`.
- Lines 16: computes `[~,s_mult]` using `[~,s_mult]=spin(sys.isotopes{1}); s_qnum=(s_mult-1)/2`.
- Lines 17: computes `inter.coupling.matrix{1,1}` using `inter.coupling.matrix{1,1}=eeqq2nqi(1.18e6,0.53,s_qnum,[0 0 0])`.
- Lines 20: computes `inter.relaxation` using `inter.relaxation={'redfield'}`.
- Lines 21: computes `inter.equilibrium` using `inter.equilibrium='zero'`.
- Lines 22: computes `inter.rlx_keep` using `inter.rlx_keep='labframe'`.
- Lines 23: computes `inter.tau_c` using `inter.tau_c={1e-9}`.
- Lines 26: computes `bas.formalism` using `bas.formalism='sphten-liouv'`.
- Lines 27: computes `bas.approximation` using `bas.approximation='none'`.
- Lines 30: computes `spin_system` using `spin_system=create(sys,inter)`.
- Lines 34: computes `R` using `R=relaxation(spin_system)`.
- Lines 37: computes `[r1,r2]` using `[r1,r2]=rlx_nqi(s_qnum,14.1*spin(sys.isotopes{1}),1.18e6,0.53,1e-9)`.
- Lines 40: computes `Lz` using `Lz=state(spin_system,'Lz',1)`.
- Lines 41: computes `Lp` using `Lp=state(spin_system,'L+',1)`.

## Implementation structure

- 14N quadrupolar relaxation in glycine in liquid state. The
- numerical output of Spinach is compared to the analytical
- equation from the textbook.
- Calculation time: seconds
- System specification
- Spin quantum number and quadrupolar tensor
- Relaxation theory
- Basis set
- Spinach housekeeping
- Relaxation superoperator
- Textbook relaxation rate expressions
- States of interest

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `spin()`, `eeqq2nqi()`, `create()`, `basis()`, `relaxation()`, `rlx_nqi()`, `state()`, `num2str()`.
