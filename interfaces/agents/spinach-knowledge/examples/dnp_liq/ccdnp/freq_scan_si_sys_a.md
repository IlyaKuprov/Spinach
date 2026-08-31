# examples/dnp_liq/ccdnp/freq_scan_si_sys_a.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/examples/dnp_liq/ccdnp/freq_scan_si_sys_a.m`
- Signature: `freq_scan_si_sys_a()`
- Total lines: 105

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
- Lines 30-31: Coordinates for anisotropic HF; implemented by `inter.coordinates={[ 0.00 0.0000 0.0000]`.
- Lines 35-36: Basis set; implemented by `bas.formalism='sphten-liouv'`.
- Lines 39-40: Disable start-up checks; implemented by `sys.disable={'hygiene'}`.
- Lines 42-43: Relaxation theory; implemented by `inter.relaxation={'redfield'}`.
- Lines 50-51: Sequence parameters; implemented by `parameters.spins={'E'}`.
- Lines 57-58: Disable excessive printing and checking; implemented by `sys.disable={'hygiene'}; sys.output='hush'`.
- Lines 60-61: Microwave frequency offset grid; implemented by `parameters.mw_frq=2*pi*linspace(-10,30,512)*1e6`.
- Lines 63-64: Magnetic field grid, Tesla; implemented by `field_grid=linspace(1,20,128)`.
- Lines 66-67: Preallocate result; implemented by `answer=zeros(numel(parameters.mw_frq),numel(field_grid),'like',1i)`.
- Lines 69-70: Loop over magnet fields; implemented by `parfor n=1:numel(field_grid)`.
- Lines 72-73: Localise parameter arrays; implemented by `localpar=parameters; localsys=sys`.
- Lines 75-76: Set the magnet field; implemented by `localsys.magnet=field_grid(n)`.
- Lines 78-79: Spinach housekeeping; implemented by `spin_system=create(localsys,inter)`.
- Lines 82-83: Relevant operators and states; implemented by `localpar.coil=state(spin_system,'Lz','1H')`.
- Lines 87-88: Isotropic thermal equilibrium; implemented by `rho_eq=equilibrium(spin_system)`.

### Control flow inferred from the code

- Line 70: `parfor` loop over `n=1:numel(field_grid)`.

### Key state/data transformations

- Lines 16: computes `sys.isotopes` using `sys.isotopes={'1H','E','E'}`.
- Lines 19: computes `inter.zeeman.eigs{1,1}` using `inter.zeeman.eigs{1,1}=[0 10 20]`.
- Lines 20: computes `inter.zeeman.euler{1,1}` using `inter.zeeman.euler{1,1}=[0 0 0]`.
- Lines 21: computes `inter.zeeman.eigs{1,2}` using `inter.zeeman.eigs{1,2}=[1.977873 1.977798 1.977792]`.
- Lines 22: computes `inter.zeeman.euler{1,2}` using `inter.zeeman.euler{1,2}=[0 0 0]`.
- Lines 23: computes `inter.zeeman.eigs{1,3}` using `inter.zeeman.eigs{1,3}=[1.977919 1.978000 1.978000]`.
- Lines 24: computes `inter.zeeman.euler{1,3}` using `inter.zeeman.euler{1,3}=[-0.59 -0.10 0.49]`.
- Lines 27: computes `inter.coupling.scalar` using `inter.coupling.scalar=cell(3,3)`.
- Lines 28: computes `inter.coupling.scalar{2,3}` using `inter.coupling.scalar{2,3}=6.2e6`.
- Lines 31: computes `inter.coordinates` using `inter.coordinates={[ 0.00 0.0000 0.0000]`.
- Lines 36: computes `bas.formalism` using `bas.formalism='sphten-liouv'`.
- Lines 37: computes `bas.approximation` using `bas.approximation='none'`.
- Lines 40: computes `sys.disable` using `sys.disable={'hygiene'}`.
- Lines 43: computes `inter.relaxation` using `inter.relaxation={'redfield'}`.
- Lines 44: computes `inter.equilibrium` using `inter.equilibrium='zero'`.
- Lines 45: computes `inter.rlx_keep` using `inter.rlx_keep='secular'`.
- Lines 46: computes `inter.temperature` using `inter.temperature=298`.
- Lines 47: computes `inter.tau_c` using `inter.tau_c={100e-12}`.

## Implementation structure

- Steady state nuclear magnetisation as a function of microwave frequency
- offset and the magnet field in a DNP experiment with two electrons con-
- nected by exchange coupling, both coupled to a nucleus by dipolar coup-
- lings. Further particulars in:
- Calculation time: seconds
- Spin system
- Zeeman interactions
- Exchange coupling
- Coordinates for anisotropic HF
- Basis set
- Disable start-up checks
- Relaxation theory

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `field_grid()`, `create()`, `basis()`, `state()`, `operator()`, `equilibrium()`, `answer()`, `liquid()`, `kfigure()`, `xlim()`, `kxlabel()`, `kylabel()`, `kcolourbar()`.
