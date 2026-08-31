# examples/dnp_liq/ccdnp/tau_scan_main_text.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/examples/dnp_liq/ccdnp/tau_scan_main_text.m`
- Signature: `tau_scan_main_text()`
- Total lines: 107

## Purpose

Steady state nuclear magnetisation as a function of microwave frequency offset and the rotational correlation time in a DNP experiment with two electrons connected by exchange coupling, both coupled to a nucleus by dipolar couplings. Further particulars in: Calculation time: seconds

## Physical / mathematical content

- Liquid-state DNP examples. The main ingredients are electron-nuclear cross-relaxation, scalar or dipolar contact mechanisms, motional spectral densities, and field/frequency dependence of polarisation transfer.
- The relaxation model is Redfield-type perturbation theory: fluctuating interactions enter through correlation functions or spectral densities and generate a linear relaxation superoperator.
- The spin physics includes through-space magnetic dipole-dipole coupling, a rank-2 anisotropic interaction with strong orientation dependence and characteristic secular/non-secular structure.

## Numerical / algorithmic content

- The implementation explicitly addresses performance engineering through parallel or GPU execution, which matters because Spinach operators can become extremely large after basis expansion or powder/spatial lifting.

## Code-derived implementation details

### Comment-guided execution stages

- Lines 15-16: Magnetic field, Tesla; implemented by `sys.magnet=14.1`.
- Lines 18-19: Spin system; implemented by `sys.isotopes={'1H','E','E'}`.
- Lines 21-22: Zeeman interactions; implemented by `inter.zeeman.eigs{1,1}=[0 10 20]`.
- Lines 29-30: Exchange coupling; implemented by `inter.coupling.scalar=cell(3,3)`.
- Lines 33-34: Cooridnates; implemented by `inter.coordinates={[ 0.000 0.000 0.000]`.
- Lines 37-38: Basis set; implemented by `bas.formalism='sphten-liouv'`.
- Lines 41-42: Turn off startup tests; implemented by `sys.disable={'hygiene'}`.
- Lines 44-45: Relaxation theory; implemented by `inter.relaxation={'redfield'}`.
- Lines 51-52: Sequence parameters; implemented by `parameters.spins={'E'}`.
- Lines 58-59: Disable excessive printing and checking; implemented by `sys.disable={'hygiene'}; sys.output='hush'`.
- Lines 61-62: Microwave frequency offset grid; implemented by `parameters.mw_frq=2*pi*linspace(-5,10,512)*1e6`.
- Lines 64-65: Correlation time grid; implemented by `tau_c=linspace(50e-12,500e-12,128)`.
- Lines 67-68: Preallocate results; implemented by `answer=zeros(numel(parameters.mw_frq),numel(tau_c),'like',1i)`.
- Lines 70-71: Loop over correlation time; implemented by `parfor n=1:numel(tau_c)`.
- Lines 73-74: Localise parameter arrays; implemented by `localpar=parameters; localint=inter`.
- Lines 76-77: Set correlation time; implemented by `localint.tau_c={tau_c(n)}`.
- Lines 79-80: Spinach housekeeping; implemented by `spin_system=create(sys,localint)`.
- Lines 83-84: Relevant operators and states; implemented by `localpar.coil=state(spin_system,'Lz','1H')`.

### Control flow inferred from the code

- Line 71: `parfor` loop over `n=1:numel(tau_c)`.

### Key state/data transformations

- Lines 16: computes `sys.magnet` using `sys.magnet=14.1`.
- Lines 19: computes `sys.isotopes` using `sys.isotopes={'1H','E','E'}`.
- Lines 22: computes `inter.zeeman.eigs{1,1}` using `inter.zeeman.eigs{1,1}=[0 10 20]`.
- Lines 23: computes `inter.zeeman.euler{1,1}` using `inter.zeeman.euler{1,1}=[0 0 0]`.
- Lines 24: computes `inter.zeeman.eigs{1,2}` using `inter.zeeman.eigs{1,2}=[2.0034 2.0038 2.0038]`.
- Lines 25: computes `inter.zeeman.euler{1,2}` using `inter.zeeman.euler{1,2}=[-0.872 -0.013 0.868]`.
- Lines 26: computes `inter.zeeman.eigs{1,3}` using `inter.zeeman.eigs{1,3}=[2.0057 2.0030 2.0030]`.
- Lines 27: computes `inter.zeeman.euler{1,3}` using `inter.zeeman.euler{1,3}=[-1.145 0.061 1.143]`.
- Lines 30: computes `inter.coupling.scalar` using `inter.coupling.scalar=cell(3,3)`.
- Lines 31: computes `inter.coupling.scalar{2,3}` using `inter.coupling.scalar{2,3}=3e6`.
- Lines 34: computes `inter.coordinates` using `inter.coordinates={[ 0.000 0.000 0.000]`.
- Lines 38: computes `bas.formalism` using `bas.formalism='sphten-liouv'`.
- Lines 39: computes `bas.approximation` using `bas.approximation='none'`.
- Lines 42: computes `sys.disable` using `sys.disable={'hygiene'}`.
- Lines 45: computes `inter.relaxation` using `inter.relaxation={'redfield'}`.
- Lines 46: computes `inter.equilibrium` using `inter.equilibrium='zero'`.
- Lines 47: computes `inter.rlx_keep` using `inter.rlx_keep='secular'`.
- Lines 48: computes `inter.temperature` using `inter.temperature=298`.

## Implementation structure

- Steady state nuclear magnetisation as a function of microwave frequency
- offset and the rotational correlation time in a DNP experiment with two
- electrons connected by exchange coupling, both coupled to a nucleus by
- dipolar couplings. Further particulars in:
- Calculation time: seconds
- Magnetic field, Tesla
- Spin system
- Zeeman interactions
- Exchange coupling
- Cooridnates
- Basis set
- Turn off startup tests

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `tau_c()`, `create()`, `basis()`, `state()`, `operator()`, `answer()`, `liquid()`, `equilibrium()`, `kfigure()`, `scale_figure()`, `kxlabel()`, `kcolourbar()`, `kylabel()`, `text()`.
