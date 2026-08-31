# examples/dnp_liq/jdnp/fig_2_tau_and_field_traject.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/examples/dnp_liq/jdnp/fig_2_tau_and_field_traject.m`
- Signature: `fig_2_tau_and_field_traject()`
- Total lines: 104

## Purpose

Time evolution plot for JDNP: proton polarisation as a function of time for specific external fields and rotational correlation times. The inter-electron exchange coupling is set to the match- ing condition at each field. Further details in: Calculation time: seconds, line-by-line plotting

## Physical / mathematical content

- Liquid-state DNP examples. The main ingredients are electron-nuclear cross-relaxation, scalar or dipolar contact mechanisms, motional spectral densities, and field/frequency dependence of polarisation transfer.

## Numerical / algorithmic content

- Time propagation is explicit. In Spinach this usually means repeated application of matrix exponentials or propagator factorizations to density operators or state vectors in Hilbert/Liouville/Fokker-Planck space.

## Code-derived implementation details

### Comment-guided execution stages

- Lines 15-16: Load the spin system; implemented by `[sys,inter,bas,parameters]=system_specification()`.
- Lines 18-19: Experiment parameters; implemented by `parameters.mw_pwr=2*pi*250e3`.
- Lines 23-24: Magnetic field grid, Tesla; implemented by `field_grid=[0.5 3.4 7.0 11.7 14.1 23.5]`.
- Lines 26-27: Correlation time grid, seconds; implemented by `tau_c=[300e-12 400e-12 500e-12 600e-12]`.
- Lines 29-30: Get a figure going; implemented by `kfigure(); scale_figure([1.75 1.0])`.
- Lines 32-33: Loop over the field grid; implemented by `for n=1:numel(field_grid)`.
- Lines 35-36: Set magnet field; implemented by `sys.magnet=field_grid(n)`.
- Lines 38-39: Set microwave offset frequency; implemented by `f_free=g2freq(parameters.g_ref,sys.magnet)`.
- Lines 43-44: Set the exchange coupling; implemented by `electron_zeeman_tensor=g2freq(inter.zeeman.matrix{2},sys.magnet)`.
- Lines 49-50: Set figure properties; implemented by `subplot(2,3,n); hold on`.
- Lines 55-56: Loop over correlation times; implemented by `for k=1:numel(tau_c)`.
- Lines 58-59: Set correlation time; implemented by `inter.tau_c={tau_c(k)}`.
- Lines 61-62: Spinach housekeeping; implemented by `spin_system=create(sys,inter)`.
- Lines 65-66: Electron operators; implemented by `Ex=operator(spin_system,'Lx','E')`.
- Lines 69-70: Isotropic thermal equilibrium; implemented by `rho_eq=equilibrium(spin_system)`.
- Lines 72-73: Hamiltonian and relaxation superoperator; implemented by `H=hamiltonian(assume(spin_system,'esr'))`.
- Lines 76-77: Add microwave terms to the Hamiltonian; implemented by `H=H+parameters.mw_pwr*Ex+parameters.mw_off*Ez`.
- Lines 79-80: Detection state: Lz on the proton; implemented by `coil=state(spin_system,'Lz','1H')`.

### Control flow inferred from the code

- Line 33: `for` loop over `n=1:numel(field_grid)`.
- Line 56: `for` loop over `k=1:numel(tau_c)`.

### Key state/data transformations

- Lines 16: computes `[sys,inter,bas,parameters]` using `[sys,inter,bas,parameters]=system_specification()`.
- Lines 19: computes `parameters.mw_pwr` using `parameters.mw_pwr=2*pi*250e3`.
- Lines 20: computes `parameters.t_step` using `parameters.t_step=1e-3`.
- Lines 21: computes `parameters.nsteps` using `parameters.nsteps=200`.
- Lines 24: computes `field_grid` using `field_grid=[0.5 3.4 7.0 11.7 14.1 23.5]`.
- Lines 27: computes `tau_c` using `tau_c=[300e-12 400e-12 500e-12 600e-12]`.
- Lines 36: computes `sys.magnet` using `sys.magnet=field_grid(n)`.
- Lines 39: computes `f_free` using `f_free=g2freq(parameters.g_ref,sys.magnet)`.
- Lines 40: computes `f_trityl` using `f_trityl=g2freq(parameters.g_trityl,sys.magnet)`.
- Lines 41: computes `parameters.mw_off` using `parameters.mw_off=2*pi*(f_trityl-f_free)`.
- Lines 44: computes `electron_zeeman_tensor` using `electron_zeeman_tensor=g2freq(inter.zeeman.matrix{2},sys.magnet)`.
- Lines 45: computes `electron_zeeman_iso` using `electron_zeeman_iso=mean(diag(electron_zeeman_tensor))`.
- Lines 46: computes `proton_zeeman_iso` using `proton_zeeman_iso=sys.magnet*spin('1H')/(2*pi)`.
- Lines 47: computes `inter.coupling.scalar{2,3}` using `inter.coupling.scalar{2,3}=electron_zeeman_iso+proton_zeeman_iso`.
- Lines 52: computes `ktitle(['$B_0$` using `ktitle(['$B_0$ = ' num2str(field_grid(n)) ' Tesla'])`.
- Lines 59: computes `inter.tau_c` using `inter.tau_c={tau_c(k)}`.
- Lines 62: computes `spin_system` using `spin_system=create(sys,inter)`.
- Lines 66: computes `Ex` using `Ex=operator(spin_system,'Lx','E')`.

## Implementation structure

- Time evolution plot for JDNP: proton polarisation as a function
- of time for specific external fields and rotational correlation
- times. The inter-electron exchange coupling is set to the match-
- ing condition at each field. Further details in:
- Calculation time: seconds, line-by-line plotting
- Load the spin system
- Experiment parameters
- Magnetic field grid, Tesla
- Correlation time grid, seconds
- Get a figure going
- Loop over the field grid
- Set magnet field

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `system_specification()`, `kfigure()`, `scale_figure()`, `field_grid()`, `g2freq()`, `spin()`, `subplot()`, `ylim()`, `ktitle()`, `num2str()`, `kylabel()`, `kxlabel()`, `tau_c()`, `create()`, `basis()`, `operator()`.
