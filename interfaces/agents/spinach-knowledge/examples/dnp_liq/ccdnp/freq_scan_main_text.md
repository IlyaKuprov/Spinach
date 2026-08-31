# examples/dnp_liq/ccdnp/freq_scan_main_text.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/examples/dnp_liq/ccdnp/freq_scan_main_text.m`
- Signature: `freq_scan_main_text()`
- Total lines: 101

## Purpose

Steady state nuclear magnetisation as a function of microwave frequency offset and the magnet field in a DNP experiment with two electrons con- nected by exchange coupling, both coupled to a nucleus by dipolar coup- lings. Further particulars in: Calculation time: seconds

## Physical / mathematical content

- Liquid-state DNP examples. The main ingredients are electron-nuclear cross-relaxation, scalar or dipolar contact mechanisms, motional spectral densities, and field/frequency dependence of polarisation transfer.
- The relaxation model is Redfield-type perturbation theory: fluctuating interactions enter through correlation functions or spectral densities and generate a linear relaxation superoperator.
- The spin physics includes through-space magnetic dipole-dipole coupling, a rank-2 anisotropic interaction with strong orientation dependence and characteristic secular/non-secular structure.

## Numerical / algorithmic content

- The implementation explicitly addresses performance engineering through parallel or GPU execution, which matters because Spinach operators can become extremely large after basis expansion or powder/spatial lifting.

## Code-derived implementation details

### Comment-guided execution stages

- Lines 15-16: Spin system; implemented by `sys.isotopes={'1H','E','E'}`.
- Lines 18-19: Zeeman interactions; implemented by `inter.zeeman.eigs{1,1}=[0 10 20]`.
- Lines 26-27: Exchange coupling; implemented by `inter.coupling.scalar=cell(3,3)`.
- Lines 30-31: Cooridnates; implemented by `inter.coordinates={[ 0.000 0.000 0.000]`.
- Lines 34-35: Basis set; implemented by `bas.formalism='sphten-liouv'`.
- Lines 38-39: Relaxation theory; implemented by `inter.relaxation={'redfield'}`.
- Lines 46-47: Sequence parameters; implemented by `parameters.spins={'E'}`.
- Lines 53-54: Disable excessive printing and checking; implemented by `sys.disable={'hygiene'}; sys.output='hush'`.
- Lines 56-57: Microwave frequency offset grid; implemented by `parameters.mw_frq=2*pi*linspace(-5,10,512)*1e6`.
- Lines 59-60: Magnetic field grid, Tesla; implemented by `field_grid=linspace(1,20,128)`.
- Lines 62-63: Preallocate result; implemented by `answer=zeros(numel(parameters.mw_frq),numel(field_grid),'like',1i)`.
- Lines 65-66: Loop over magnet fields; implemented by `parfor n=1:numel(field_grid)`.
- Lines 68-69: Localise parameter arrays; implemented by `localpar=parameters; localsys=sys`.
- Lines 71-72: Set the magnet field; implemented by `localsys.magnet=field_grid(n)`.
- Lines 74-75: Spinach housekeeping; implemented by `spin_system=create(localsys,inter)`.
- Lines 78-79: Relevant operators and states; implemented by `localpar.coil=state(spin_system,'Lz','1H')`.
- Lines 83-84: Isotropic thermal equilibrium; implemented by `rho_eq=equilibrium(spin_system)`.
- Lines 86-87: Steady state simulation; implemented by `answer(:,n)=liquid(spin_system,@dnp_freq_scan,localpar,'esr')`.

### Control flow inferred from the code

- Line 66: `parfor` loop over `n=1:numel(field_grid)`.

### Key state/data transformations

- Lines 16: computes `sys.isotopes` using `sys.isotopes={'1H','E','E'}`.
- Lines 19: computes `inter.zeeman.eigs{1,1}` using `inter.zeeman.eigs{1,1}=[0 10 20]`.
- Lines 20: computes `inter.zeeman.euler{1,1}` using `inter.zeeman.euler{1,1}=[0 0 0]`.
- Lines 21: computes `inter.zeeman.eigs{1,2}` using `inter.zeeman.eigs{1,2}=[2.0034 2.0038 2.0038]`.
- Lines 22: computes `inter.zeeman.euler{1,2}` using `inter.zeeman.euler{1,2}=[-0.872 -0.013 0.868]`.
- Lines 23: computes `inter.zeeman.eigs{1,3}` using `inter.zeeman.eigs{1,3}=[2.0057 2.0030 2.0030]`.
- Lines 24: computes `inter.zeeman.euler{1,3}` using `inter.zeeman.euler{1,3}=[-1.145 0.061 1.143]`.
- Lines 27: computes `inter.coupling.scalar` using `inter.coupling.scalar=cell(3,3)`.
- Lines 28: computes `inter.coupling.scalar{2,3}` using `inter.coupling.scalar{2,3}=3e6`.
- Lines 31: computes `inter.coordinates` using `inter.coordinates={[ 0.000 0.000 0.000]`.
- Lines 35: computes `bas.formalism` using `bas.formalism='sphten-liouv'`.
- Lines 36: computes `bas.approximation` using `bas.approximation='none'`.
- Lines 39: computes `inter.relaxation` using `inter.relaxation={'redfield'}`.
- Lines 40: computes `inter.equilibrium` using `inter.equilibrium='zero'`.
- Lines 41: computes `inter.rlx_keep` using `inter.rlx_keep='secular'`.
- Lines 42: computes `inter.temperature` using `inter.temperature=298`.
- Lines 43: computes `inter.tau_c` using `inter.tau_c={100e-12}`.
- Lines 44: computes `sys.tols.rlx_integration` using `sys.tols.rlx_integration=1e-10`.

## Implementation structure

- Steady state nuclear magnetisation as a function of microwave frequency
- offset and the magnet field in a DNP experiment with two electrons con-
- nected by exchange coupling, both coupled to a nucleus by dipolar coup-
- lings. Further particulars in:
- Calculation time: seconds
- Spin system
- Zeeman interactions
- Exchange coupling
- Cooridnates
- Basis set
- Relaxation theory
- Sequence parameters

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `field_grid()`, `create()`, `basis()`, `state()`, `operator()`, `equilibrium()`, `answer()`, `liquid()`, `kfigure()`, `xlim()`, `kxlabel()`, `kylabel()`, `kcolourbar()`.
