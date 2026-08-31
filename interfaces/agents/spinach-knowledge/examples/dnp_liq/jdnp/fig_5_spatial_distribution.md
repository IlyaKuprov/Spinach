# examples/dnp_liq/jdnp/fig_5_spatial_distribution.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/examples/dnp_liq/jdnp/fig_5_spatial_distribution.m`
- Signature: `fig_5_spatial_distribution()`
- Total lines: 156

## Purpose

An illustration of the fact that JDNP effect does not vanish on position and orientation averaging in liquid phase. The simula- tion shows proton polarisation at 20 ms for different proton lo- cations around a radical pair with inter-electron exchange cou- pling chosen to acheive the JDNP effect. Details in: Calculation time: seconds

## Physical / mathematical content

- Liquid-state DNP examples. The main ingredients are electron-nuclear cross-relaxation, scalar or dipolar contact mechanisms, motional spectral densities, and field/frequency dependence of polarisation transfer.

## Numerical / algorithmic content

- Time propagation is explicit. In Spinach this usually means repeated application of matrix exponentials or propagator factorizations to density operators or state vectors in Hilbert/Liouville/Fokker-Planck space.
- The implementation explicitly addresses performance engineering through parallel or GPU execution, which matters because Spinach operators can become extremely large after basis expansion or powder/spatial lifting.

## Code-derived implementation details

### Comment-guided execution stages

- Lines 16-17: Load the spin system; implemented by `[sys,inter,bas,parameters]=system_specification()`.
- Lines 19-20: Set magnet field; implemented by `sys.magnet=14.08`.
- Lines 22-23: Experiment parameters; implemented by `parameters.mw_pwr=2*pi*250e3`.
- Lines 26-27: Set microwave offset frequency; implemented by `f_free=g2freq(parameters.g_ref,sys.magnet)`.
- Lines 31-32: Set the exchange coupling; implemented by `electron_zeeman_tensor=g2freq(inter.zeeman.matrix{2},sys.magnet)`.
- Lines 37-38: Specify coordinate arrays; implemented by `X=linspace(-30,30,30); Y=linspace(-30,30,30); Z=linspace(-30,30,30)`.
- Lines 40-41: Get a figure going; implemented by `kfigure(); scale_figure([1.75 1.0])`.
- Lines 43-44: Preallocate answer arrays; implemented by `xy_dnp=nan(numel(X),numel(Y))`.
- Lines 47-48: XY plane scan, Z=0; implemented by `for n=1:numel(X)`.
- Lines 50-51: Parallel inner loop; implemented by `parfor k=1:numel(Y)`.
- Lines 53-54: Localise arrays; implemented by `localpar=parameters; localint=inter`.
- Lines 56-57: Set proton position; implemented by `localint.coordinates{1,1}=[X(n) Y(k) 0.0]`.
- Lines 59-60: Spinach housekeeping; implemented by `spin_system=create(sys,localint)`.
- Lines 63-64: Electron operators; implemented by `Ex=operator(spin_system,'Lx','E')`.
- Lines 67-68: Detection state; implemented by `Nz=state(spin_system,'Lz','1H')`.
- Lines 70-71: Isotropic thermal equilibrium; implemented by `rho_eq=equilibrium(spin_system)`.
- Lines 73-74: Hamiltonian and relaxation superoperator; implemented by `H=hamiltonian(assume(spin_system,'esr'))`.
- Lines 77-78: Add microwave terms to the Hamiltonian; implemented by `H=H+parameters.mw_pwr*Ex+parameters.mw_off*Ez`.

### Control flow inferred from the code

- Line 48: `for` loop over `n=1:numel(X)`.
- Line 51: `parfor` loop over `k=1:numel(Y)`.
- Line 102: `for` loop over `n=1:numel(X)`.
- Line 105: `parfor` loop over `k=1:numel(Z)`.

### Key state/data transformations

- Lines 17: computes `[sys,inter,bas,parameters]` using `[sys,inter,bas,parameters]=system_specification()`.
- Lines 20: computes `sys.magnet` using `sys.magnet=14.08`.
- Lines 23: computes `parameters.mw_pwr` using `parameters.mw_pwr=2*pi*250e3`.
- Lines 24: computes `parameters.t_evol` using `parameters.t_evol=20e-3`.
- Lines 27: computes `f_free` using `f_free=g2freq(parameters.g_ref,sys.magnet)`.
- Lines 28: computes `f_trityl` using `f_trityl=g2freq(parameters.g_trityl,sys.magnet)`.
- Lines 29: computes `parameters.mw_off` using `parameters.mw_off=2*pi*(f_trityl-f_free)`.
- Lines 32: computes `electron_zeeman_tensor` using `electron_zeeman_tensor=g2freq(inter.zeeman.matrix{2},sys.magnet)`.
- Lines 33: computes `electron_zeeman_iso` using `electron_zeeman_iso=mean(diag(electron_zeeman_tensor))`.
- Lines 34: computes `proton_zeeman_iso` using `proton_zeeman_iso=sys.magnet*spin('1H')/(2*pi)`.
- Lines 35: computes `inter.coupling.scalar{2,3}` using `inter.coupling.scalar{2,3}=electron_zeeman_iso+proton_zeeman_iso`.
- Lines 38: computes `X` using `X=linspace(-30,30,30); Y=linspace(-30,30,30); Z=linspace(-30,30,30)`.
- Lines 44: computes `xy_dnp` using `xy_dnp=nan(numel(X),numel(Y))`.
- Lines 45: computes `xz_dnp` using `xz_dnp=nan(numel(X),numel(Z))`.
- Lines 54: computes `localpar` using `localpar=parameters; localint=inter`.
- Lines 57: computes `localint.coordinates{1,1}` using `localint.coordinates{1,1}=[X(n) Y(k) 0.0]`.
- Lines 60: computes `spin_system` using `spin_system=create(sys,localint)`.
- Lines 64: computes `Ex` using `Ex=operator(spin_system,'Lx','E')`.

## Implementation structure

- An illustration of the fact that JDNP effect does not vanish on
- position and orientation averaging in liquid phase. The simula-
- tion shows proton polarisation at 20 ms for different proton lo-
- cations around a radical pair with inter-electron exchange cou-
- pling chosen to acheive the JDNP effect. Details in:
- Calculation time: seconds
- Load the spin system
- Set magnet field
- Experiment parameters
- Set microwave offset frequency
- Set the exchange coupling
- Specify coordinate arrays

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `system_specification()`, `g2freq()`, `spin()`, `kfigure()`, `scale_figure()`, `nan()`, `create()`, `basis()`, `operator()`, `state()`, `equilibrium()`, `hamiltonian()`, `assume()`, `relaxation()`, `evolution()`, `xy_dnp()`.
