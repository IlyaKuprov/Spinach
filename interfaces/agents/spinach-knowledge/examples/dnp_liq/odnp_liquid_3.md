# examples/dnp_liq/odnp_liquid_3.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/examples/dnp_liq/odnp_liquid_3.m`
- Signature: `odnp_liquid_3()`
- Total lines: 95

## Purpose

Steady state nuclear magnetisation as a function of microwave frequency offset and the magnet field in a DNP experiment with an electron and a nucleus connected by a hyperfine coupling. A g-hyperfine cross-correla- tion effect is visible at high field. The steady state is computed by setting the time derivative to zero in the inhomogeneous master equation, and solving the resulting algebraic equation for the steady s

## Physical / mathematical content

- Liquid-state DNP examples. The main ingredients are electron-nuclear cross-relaxation, scalar or dipolar contact mechanisms, motional spectral densities, and field/frequency dependence of polarisation transfer.
- The relaxation model is Redfield-type perturbation theory: fluctuating interactions enter through correlation functions or spectral densities and generate a linear relaxation superoperator.
- The spin physics includes through-space magnetic dipole-dipole coupling, a rank-2 anisotropic interaction with strong orientation dependence and characteristic secular/non-secular structure.

## Numerical / algorithmic content

- The implementation explicitly addresses performance engineering through parallel or GPU execution, which matters because Spinach operators can become extremely large after basis expansion or powder/spatial lifting.

## Code-derived implementation details

### Comment-guided execution stages

- Lines 17-18: Spin system; implemented by `sys.isotopes={'1H','E'}`.
- Lines 20-21: Anisotropic Zeeman interactions; implemented by `inter.zeeman.eigs{1}=[15 5 -20]`.
- Lines 26-27: Isotropic hyperfine coupling; implemented by `inter.coupling.scalar{1,2}=20e6`.
- Lines 30-31: Coordinates for dipolar coupling; implemented by `inter.coordinates={[0.0 0.0 0.0]`.
- Lines 34-35: Basis set; implemented by `bas.formalism='sphten-liouv'`.
- Lines 38-39: Relaxation theory; implemented by `inter.relaxation={'redfield'}`.
- Lines 46-47: Sequence parameters; implemented by `parameters.spins={'E'}`.
- Lines 53-54: Turn off excessive output; implemented by `sys.disable={'hygiene'}; sys.output='hush'`.
- Lines 56-57: Microwave frequency offset grid; implemented by `parameters.mw_frq=2*pi*linspace(-15,15,512)*1e6`.
- Lines 59-60: Magnetic field grid, Tesla; implemented by `field_grid=linspace(1,10,64)`.
- Lines 62-63: Preallocate result; implemented by `answer=zeros(numel(parameters.mw_frq),numel(field_grid))`.
- Lines 65-66: Loop over magnet fields; implemented by `parfor n=1:numel(field_grid)`.
- Lines 68-69: Localise parameter arrays; implemented by `localsys=sys; locpar=parameters`.
- Lines 71-72: Set the magnet field; implemented by `localsys.magnet=field_grid(n)`.
- Lines 74-75: Spinach housekeeping; implemented by `spin_system=create(localsys,inter)`.
- Lines 78-79: Relevant operators and states; implemented by `locpar.coil=state(spin_system,'Lz','1H')`.
- Lines 83-84: Steady state simulation; implemented by `answer(:,n)=liquid(spin_system,@dnp_freq_scan,locpar,'esr')`.
- Lines 88-89: Plotting; implemented by `kfigure(); imagesc(field_grid,(parameters.mw_frq/(2*pi*1e6)),real(answer))`.

### Control flow inferred from the code

- Line 66: `parfor` loop over `n=1:numel(field_grid)`.

### Key state/data transformations

- Lines 18: computes `sys.isotopes` using `sys.isotopes={'1H','E'}`.
- Lines 21: computes `inter.zeeman.eigs{1}` using `inter.zeeman.eigs{1}=[15 5 -20]`.
- Lines 22: computes `inter.zeeman.eigs{2}` using `inter.zeeman.eigs{2}=[2.00210 2.00250 2.00290]`.
- Lines 23: computes `inter.zeeman.euler{1}` using `inter.zeeman.euler{1}=[0.00 0.00 0.00]`.
- Lines 24: computes `inter.zeeman.euler{2}` using `inter.zeeman.euler{2}=[pi/3 pi/4 pi/5]`.
- Lines 27: computes `inter.coupling.scalar{1,2}` using `inter.coupling.scalar{1,2}=20e6`.
- Lines 28: computes `inter.coupling.scalar{2,2}` using `inter.coupling.scalar{2,2}=0`.
- Lines 31: computes `inter.coordinates` using `inter.coordinates={[0.0 0.0 0.0]`.
- Lines 35: computes `bas.formalism` using `bas.formalism='sphten-liouv'`.
- Lines 36: computes `bas.approximation` using `bas.approximation='none'`.
- Lines 39: computes `inter.relaxation` using `inter.relaxation={'redfield'}`.
- Lines 40: computes `inter.equilibrium` using `inter.equilibrium='zero'`.
- Lines 41: computes `inter.rlx_keep` using `inter.rlx_keep='secular'`.
- Lines 42: computes `inter.temperature` using `inter.temperature=298`.
- Lines 43: computes `inter.tau_c` using `inter.tau_c={10e-12}`.
- Lines 44: computes `sys.tols.rlx_integration` using `sys.tols.rlx_integration=1e-10`.
- Lines 47: computes `parameters.spins` using `parameters.spins={'E'}`.
- Lines 48: computes `parameters.mw_pwr` using `parameters.mw_pwr=2*pi*500e3`.

## Implementation structure

- Steady state nuclear magnetisation as a function of microwave frequency
- offset and the magnet field in a DNP experiment with an electron and a
- nucleus connected by a hyperfine coupling. A g-hyperfine cross-correla-
- tion effect is visible at high field.
- The steady state is computed by setting the time derivative to zero in
- the inhomogeneous master equation, and solving the resulting algebraic
- equation for the steady state density matrix.
- Calculation time: minutes.
- Spin system
- Anisotropic Zeeman interactions
- Isotropic hyperfine coupling
- Coordinates for dipolar coupling

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `field_grid()`, `create()`, `basis()`, `state()`, `operator()`, `answer()`, `liquid()`, `kfigure()`, `xlim()`, `kxlabel()`, `kylabel()`, `kcolourbar()`.
