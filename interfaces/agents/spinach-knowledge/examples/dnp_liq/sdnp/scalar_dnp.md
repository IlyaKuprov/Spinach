# examples/dnp_liq/sdnp/scalar_dnp.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/examples/dnp_liq/sdnp/scalar_dnp.m`
- Signature: `scalar_dnp()`
- Total lines: 116

## Purpose

Field dependence of the couping factor between 13C of CHCl3 and the electron spin of a nitroxide radical. Further particulars here: Experimental data from: Calculation time: seconds

## Physical / mathematical content

- Liquid-state DNP examples. The main ingredients are electron-nuclear cross-relaxation, scalar or dipolar contact mechanisms, motional spectral densities, and field/frequency dependence of polarisation transfer.
- The relaxation model is Redfield-type perturbation theory: fluctuating interactions enter through correlation functions or spectral densities and generate a linear relaxation superoperator.
- The spin physics includes through-space magnetic dipole-dipole coupling, a rank-2 anisotropic interaction with strong orientation dependence and characteristic secular/non-secular structure.
- Chemical-shift anisotropy is present: shielding is treated as a second-rank tensor whose orientation relative to the field or rotor axis modulates line shapes and transfer dynamics.

## Numerical / algorithmic content

- The implementation explicitly addresses performance engineering through parallel or GPU execution, which matters because Spinach operators can become extremely large after basis expansion or powder/spatial lifting.

## Code-derived implementation details

### Comment-guided execution stages

- Lines 17-18: Spin system; implemented by `sys.isotopes={'E','13C'}`.
- Lines 20-21: Coordinates for dipolar Redfield; implemented by `inter.coordinates{1}=[0.0 0.0 0.0]`.
- Lines 24-25: Zeeman interactions for CSA/g-aniso Redfield; implemented by `inter.zeeman.eigs{1}=[2.0029 2.0065 2.0098]`.
- Lines 30-31: Static hyperfine coupling; implemented by `inter.coupling.scalar=cell(2,2)`.
- Lines 34-35: Formalism and approximation; implemented by `bas.formalism='sphten-liouv'`.
- Lines 38-39: Relaxation theories; implemented by `inter.relaxation={'SRFK','redfield','t1_t2'}`.
- Lines 41-42: Electron R1 and R2 for empirical T1/T2; implemented by `inter.r1_rates{1}=2e6`.
- Lines 45-46: Nuclear R1 and R2 for empirical T1/T2; implemented by `inter.r1_rates{2}=0.25`.
- Lines 49-50: Correlation time for rotational Redfield; implemented by `inter.tau_c={30e-12}`.
- Lines 52-53: Parameters of scalar collisional Redfield; implemented by `inter.srfk_tau_c={[0.62, 30.0e-12],`.
- Lines 58-59: Temperature and thermalisation; implemented by `inter.temperature=298`.
- Lines 62-63: Secular approximation; implemented by `inter.rlx_keep='secular'`.
- Lines 65-66: Magnetic field axis; implemented by `b_vector=[logspace(-2,1,30) 15 20 25 30]`.
- Lines 68-69: Prevent excessive output; implemented by `sys.disable={'hygiene'}; sys.output='hush'`.
- Lines 71-72: Loop over magnet fields; implemented by `Rx=zeros(1,34); R1n=zeros(1,34)`.
- Lines 75-76: Localise and set the magnet field; implemented by `locsys=sys; locsys.magnet=b_vector(n)`.
- Lines 78-79: Run Spinach housekeeping; implemented by `spin_system = create(locsys,inter)`.
- Lines 82-83: Get the relaxation superoperator; implemented by `R=relaxation(spin_system)`.

### Control flow inferred from the code

- Line 73: `parfor` loop over `n=1:numel(b_vector)`.

### Key state/data transformations

- Lines 18: computes `sys.isotopes` using `sys.isotopes={'E','13C'}`.
- Lines 21: computes `inter.coordinates{1}` using `inter.coordinates{1}=[0.0 0.0 0.0]`.
- Lines 22: computes `inter.coordinates{2}` using `inter.coordinates{2}=[0.0 0.0 3.1]`.
- Lines 25: computes `inter.zeeman.eigs{1}` using `inter.zeeman.eigs{1}=[2.0029 2.0065 2.0098]`.
- Lines 26: computes `inter.zeeman.eigs{2}` using `inter.zeeman.eigs{2}=[100 120 120]`.
- Lines 27: computes `inter.zeeman.euler{1}` using `inter.zeeman.euler{1}=[0 0 0]`.
- Lines 28: computes `inter.zeeman.euler{2}` using `inter.zeeman.euler{2}=[0 0 0]`.
- Lines 31: computes `inter.coupling.scalar` using `inter.coupling.scalar=cell(2,2)`.
- Lines 32: computes `inter.coupling.scalar{1,2}` using `inter.coupling.scalar{1,2}=2e6`.
- Lines 35: computes `bas.formalism` using `bas.formalism='sphten-liouv'`.
- Lines 36: computes `bas.approximation` using `bas.approximation='none'`.
- Lines 39: computes `inter.relaxation` using `inter.relaxation={'SRFK','redfield','t1_t2'}`.
- Lines 42: computes `inter.r1_rates{1}` using `inter.r1_rates{1}=2e6`.
- Lines 43: computes `inter.r2_rates{1}` using `inter.r2_rates{1}=5e6`.
- Lines 46: computes `inter.r1_rates{2}` using `inter.r1_rates{2}=0.25`.
- Lines 47: computes `inter.r2_rates{2}` using `inter.r2_rates{2}=3.00`.
- Lines 50: computes `inter.tau_c` using `inter.tau_c={30e-12}`.
- Lines 53: computes `inter.srfk_tau_c` using `inter.srfk_tau_c={[0.62, 30.0e-12],`.

## Implementation structure

- Field dependence of the couping factor between 13C of CHCl3 and the
- electron spin of a nitroxide radical. Further particulars here:
- Experimental data from:
- Calculation time: seconds
- Spin system
- Coordinates for dipolar Redfield
- Zeeman interactions for CSA/g-aniso Redfield
- Static hyperfine coupling
- Formalism and approximation
- Relaxation theories
- Electron R1 and R2 for empirical T1/T2
- Nuclear R1 and R2 for empirical T1/T2

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `b_vector()`, `create()`, `basis()`, `relaxation()`, `state()`, `R1n()`, `kfigure()`, `expt_data()`, `set()`, `kxlabel()`, `kylabel()`.
