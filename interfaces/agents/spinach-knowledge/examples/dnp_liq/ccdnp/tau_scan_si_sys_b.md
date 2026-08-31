# examples/dnp_liq/ccdnp/tau_scan_si_sys_b.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/examples/dnp_liq/ccdnp/tau_scan_si_sys_b.m`
- Signature: `tau_scan_si_sys_b()`
- Total lines: 104

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
- Lines 33-34: coordinates; implemented by `inter.coordinates={[ 0.000 0.000 0.000]`.
- Lines 38-39: Basis set; implemented by `bas.formalism='sphten-liouv'`.
- Lines 42-43: Relaxation theory; implemented by `inter.relaxation={'redfield'}`.
- Lines 49-50: Sequence parameters; implemented by `parameters.spins={'E'}`.
- Lines 56-57: Disable excessive printing and checking; implemented by `sys.disable={'hygiene'}; sys.output='hush'`.
- Lines 59-60: Microwave frequency offset grid; implemented by `parameters.mw_frq=2*pi*linspace(-15,15,512)*1e6`.
- Lines 62-63: Correlation time grid; implemented by `tau_c=linspace(50e-12,500e-12,128)`.
- Lines 65-66: Preallocate results; implemented by `answer=zeros(numel(parameters.mw_frq),numel(tau_c),'like',1i)`.
- Lines 68-69: Loop over correlation time; implemented by `parfor n=1:numel(tau_c)`.
- Lines 71-72: Localise parameter arrays; implemented by `localpar=parameters; localint=inter`.
- Lines 74-75: Set correlation time; implemented by `localint.tau_c={tau_c(n)}`.
- Lines 77-78: Spinach housekeeping; implemented by `spin_system=create(sys,localint)`.
- Lines 81-82: Relevant operators and states; implemented by `localpar.coil=state(spin_system,'Lz','1H')`.
- Lines 86-87: Steady state simulation; implemented by `answer(:,n)=liquid(spin_system,@dnp_freq_scan,localpar,'esr')`.

### Control flow inferred from the code

- Line 69: `parfor` loop over `n=1:numel(tau_c)`.

### Key state/data transformations

- Lines 16: computes `sys.magnet` using `sys.magnet=14.1`.
- Lines 19: computes `sys.isotopes` using `sys.isotopes={'1H','E','E'}`.
- Lines 22: computes `inter.zeeman.eigs{1,1}` using `inter.zeeman.eigs{1,1}=[0 10 20]`.
- Lines 23: computes `inter.zeeman.euler{1,1}` using `inter.zeeman.euler{1,1}=[0 0 0]`.
- Lines 24: computes `inter.zeeman.eigs{1,2}` using `inter.zeeman.eigs{1,2}=[1.9778 1.9776 1.9776]`.
- Lines 25: computes `inter.zeeman.euler{1,2}` using `inter.zeeman.euler{1,2}=[-0.180 0.017 0.194]`.
- Lines 26: computes `inter.zeeman.eigs{1,3}` using `inter.zeeman.eigs{1,3}=[2.0068 2.0038 2.0038]`.
- Lines 27: computes `inter.zeeman.euler{1,3}` using `inter.zeeman.euler{1,3}=[0.632 0.783 1.086]`.
- Lines 30: computes `inter.coupling.scalar` using `inter.coupling.scalar=cell(3,3)`.
- Lines 31: computes `inter.coupling.scalar{2,3}` using `inter.coupling.scalar{2,3}=5e6`.
- Lines 34: computes `inter.coordinates` using `inter.coordinates={[ 0.000 0.000 0.000]`.
- Lines 39: computes `bas.formalism` using `bas.formalism='sphten-liouv'`.
- Lines 40: computes `bas.approximation` using `bas.approximation='none'`.
- Lines 43: computes `inter.relaxation` using `inter.relaxation={'redfield'}`.
- Lines 44: computes `inter.equilibrium` using `inter.equilibrium='zero'`.
- Lines 45: computes `inter.rlx_keep` using `inter.rlx_keep='secular'`.
- Lines 46: computes `inter.temperature` using `inter.temperature=298`.
- Lines 47: computes `sys.tols.rlx_integration` using `sys.tols.rlx_integration=1e-10`.

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
- coordinates
- Basis set
- Relaxation theory

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `tau_c()`, `create()`, `basis()`, `state()`, `operator()`, `answer()`, `liquid()`, `equilibrium()`, `kfigure()`, `scale_figure()`, `kxlabel()`, `kcolourbar()`, `kylabel()`.
