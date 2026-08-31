# examples/dnp_liq/jdnp/system_specification.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/examples/dnp_liq/jdnp/system_specification.m`
- Signature: `[sys,inter,bas,parameters]=system_specification()`
- Total lines: 50

## Purpose

Parameters of the 2e1n system used for the simulations reported in https://doi.org/10.1039/d1cp04186j

## Physical / mathematical content

- Liquid-state DNP examples. The main ingredients are electron-nuclear cross-relaxation, scalar or dipolar contact mechanisms, motional spectral densities, and field/frequency dependence of polarisation transfer.
- The relaxation model is Redfield-type perturbation theory: fluctuating interactions enter through correlation functions or spectral densities and generate a linear relaxation superoperator.

## Numerical / algorithmic content


## Code-derived implementation details

### Comment-guided execution stages

- Lines 9-10: Spin system; implemented by `sys.isotopes={'1H','E','E'}`.
- Lines 12-13: Nuclear chemical shift tensor; implemented by `inter.zeeman.matrix{1}=diag([5 10 20])`.
- Lines 15-16: Electron g-tensors -axial along Z; implemented by `inter.zeeman.matrix{2}=diag([2.0032 2.0032 2.0026])`.
- Lines 19-20: Set coordinates; implemented by `inter.coordinates{1,1}=[-3.00 0.50 1.30]`.
- Lines 24-25: Empty scalar coupling array; implemented by `inter.coupling.scalar=cell(3,3)`.
- Lines 27-28: Relaxation theory; implemented by `inter.relaxation={'redfield','SRFK'}`.
- Lines 37-38: Basis set; implemented by `bas.formalism='sphten-liouv'`.
- Lines 41-42: Tolerance settings; implemented by `sys.tols.rlx_integration=1e-10`.
- Lines 46-47: Reference g-factors; implemented by `parameters.g_ref=2.00231930436256`.

### Key state/data transformations

- Lines 10: computes `sys.isotopes` using `sys.isotopes={'1H','E','E'}`.
- Lines 13: computes `inter.zeeman.matrix{1}` using `inter.zeeman.matrix{1}=diag([5 10 20])`.
- Lines 16: computes `inter.zeeman.matrix{2}` using `inter.zeeman.matrix{2}=diag([2.0032 2.0032 2.0026])`.
- Lines 17: computes `inter.zeeman.matrix{3}` using `inter.zeeman.matrix{3}=diag([2.0032 2.0032 2.0026])`.
- Lines 20: computes `inter.coordinates{1,1}` using `inter.coordinates{1,1}=[-3.00 0.50 1.30]`.
- Lines 21: computes `inter.coordinates{2,1}` using `inter.coordinates{2,1}=[ 0.00 0.00 -9.37]`.
- Lines 22: computes `inter.coordinates{3,1}` using `inter.coordinates{3,1}=[ 0.00 0.00 +9.37]`.
- Lines 25: computes `inter.coupling.scalar` using `inter.coupling.scalar=cell(3,3)`.
- Lines 28: computes `inter.relaxation` using `inter.relaxation={'redfield','SRFK'}`.
- Lines 29: computes `inter.equilibrium` using `inter.equilibrium='dibari'`.
- Lines 30: computes `inter.rlx_keep` using `inter.rlx_keep='labframe'`.
- Lines 31: computes `inter.temperature` using `inter.temperature=298`.
- Lines 32: computes `inter.tau_c` using `inter.tau_c={500e-12}`.
- Lines 33: computes `inter.srfk_tau_c` using `inter.srfk_tau_c={[1.0 1e-12]}`.
- Lines 34: computes `inter.srfk_mdepth` using `inter.srfk_mdepth=cell(3)`.
- Lines 35: computes `inter.srfk_mdepth{2,3}` using `inter.srfk_mdepth{2,3}=3e9`.
- Lines 38: computes `bas.formalism` using `bas.formalism='sphten-liouv'`.
- Lines 39: computes `bas.approximation` using `bas.approximation='none'`.

## Implementation structure

- Parameters of the 2e1n system used for the simulations reported
- in https://doi.org/10.1039/d1cp04186j
- Spin system
- Nuclear chemical shift tensor
- Electron g-tensors -axial along Z
- Set coordinates
- Empty scalar coupling array
- Relaxation theory
- Basis set
- Tolerance settings
- Reference g-factors
