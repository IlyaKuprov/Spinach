# examples/dnp_liq/jdnp/fig_3_time_dep_bot_row.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/examples/dnp_liq/jdnp/fig_3_time_dep_bot_row.m`
- Signature: `fig_3_time_dep_bot_row()`
- Total lines: 86

## Purpose

Time evolution plot for JDNP: proton polarisation as a function of time for specific external fields. The inter-electron exchan- ge coupling is set to the matching condition at each field. See also the "top row" simulation where one of the electrons is re- moved to demostrate that JDNP vanishes. Further details in: Calculation time: seconds, line-by-line plotting

## Physical / mathematical content

- Liquid-state DNP examples. The main ingredients are electron-nuclear cross-relaxation, scalar or dipolar contact mechanisms, motional spectral densities, and field/frequency dependence of polarisation transfer.

## Numerical / algorithmic content

- Time propagation is explicit. In Spinach this usually means repeated application of matrix exponentials or propagator factorizations to density operators or state vectors in Hilbert/Liouville/Fokker-Planck space.

## Code-derived implementation details

### Comment-guided execution stages

- Lines 16-17: Load the spin system; implemented by `[sys,inter,bas,parameters]=system_specification()`.
- Lines 19-20: Experiment parameters; implemented by `parameters.mw_pwr=2*pi*250e3`.
- Lines 24-25: Magnetic field grid, Tesla; implemented by `field_grid=[0.034 0.34 3.4]`.
- Lines 27-28: Get a figure going; implemented by `kfigure(); scale_figure([1.75 0.5])`.
- Lines 30-31: Loop over the fields; implemented by `for n=1:numel(field_grid)`.
- Lines 33-34: Set magnet field; implemented by `sys.magnet=field_grid(n)`.
- Lines 36-37: Set microwave offset frequency; implemented by `f_free=g2freq(parameters.g_ref,sys.magnet)`.
- Lines 41-42: Set the exchange coupling; implemented by `electron_zeeman_tensor=g2freq(inter.zeeman.matrix{2},sys.magnet)`.
- Lines 47-48: Spinach housekeeping; implemented by `spin_system=create(sys,inter)`.
- Lines 51-52: Electron operators; implemented by `Ex=operator(spin_system,'Lx','E')`.
- Lines 55-56: Isotropic thermal equilibrium; implemented by `rho_eq=equilibrium(spin_system)`.
- Lines 58-59: Hamiltonian and relaxation superoperator; implemented by `H=hamiltonian(assume(spin_system,'esr'))`.
- Lines 62-63: Add microwave terms to the Hamiltonian; implemented by `H=H+parameters.mw_pwr*Ex+parameters.mw_off*Ez`.
- Lines 65-66: Detection state: Lz on the proton; implemented by `coil=state(spin_system,'Lz','1H')`.
- Lines 68-70: Run the time evolution and normalise to thermal equilibrium; implemented by `answer=evolution(spin_system,H+1i*R,coil,rho_eq,parameters.t_step, parameters.nsteps,'observable')`.
- Lines 73-75: Time axis generation; implemented by `t_axis=linspace(0,parameters.t_step*parameters.nsteps, parameters.nsteps+1)`.
- Lines 77-78: Plotting; implemented by `subplot(1,3,n); plot(t_axis,answer); kgrid; box on`.

### Control flow inferred from the code

- Line 31: `for` loop over `n=1:numel(field_grid)`.

### Key state/data transformations

- Lines 17: computes `[sys,inter,bas,parameters]` using `[sys,inter,bas,parameters]=system_specification()`.
- Lines 20: computes `parameters.mw_pwr` using `parameters.mw_pwr=2*pi*250e3`.
- Lines 21: computes `parameters.t_step` using `parameters.t_step=1e-3`.
- Lines 22: computes `parameters.nsteps` using `parameters.nsteps=300`.
- Lines 25: computes `field_grid` using `field_grid=[0.034 0.34 3.4]`.
- Lines 34: computes `sys.magnet` using `sys.magnet=field_grid(n)`.
- Lines 37: computes `f_free` using `f_free=g2freq(parameters.g_ref,sys.magnet)`.
- Lines 38: computes `f_trityl` using `f_trityl=g2freq(parameters.g_trityl,sys.magnet)`.
- Lines 39: computes `parameters.mw_off` using `parameters.mw_off=2*pi*(f_trityl-f_free)`.
- Lines 42: computes `electron_zeeman_tensor` using `electron_zeeman_tensor=g2freq(inter.zeeman.matrix{2},sys.magnet)`.
- Lines 43: computes `electron_zeeman_iso` using `electron_zeeman_iso=mean(diag(electron_zeeman_tensor))`.
- Lines 44: computes `proton_zeeman_iso` using `proton_zeeman_iso=sys.magnet*spin('1H')/(2*pi)`.
- Lines 45: computes `inter.coupling.scalar{2,3}` using `inter.coupling.scalar{2,3}=electron_zeeman_iso+proton_zeeman_iso`.
- Lines 48: computes `spin_system` using `spin_system=create(sys,inter)`.
- Lines 52: computes `Ex` using `Ex=operator(spin_system,'Lx','E')`.
- Lines 53: computes `Ez` using `Ez=operator(spin_system,'Lz','E')`.
- Lines 56: computes `rho_eq` using `rho_eq=equilibrium(spin_system)`.
- Lines 59: computes `H` using `H=hamiltonian(assume(spin_system,'esr'))`.

## Implementation structure

- Time evolution plot for JDNP: proton polarisation as a function
- of time for specific external fields. The inter-electron exchan-
- ge coupling is set to the matching condition at each field. See
- also the "top row" simulation where one of the electrons is re-
- moved to demostrate that JDNP vanishes. Further details in:
- Calculation time: seconds, line-by-line plotting
- Load the spin system
- Experiment parameters
- Magnetic field grid, Tesla
- Get a figure going
- Loop over the fields
- Set magnet field

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `system_specification()`, `kfigure()`, `scale_figure()`, `field_grid()`, `g2freq()`, `spin()`, `create()`, `basis()`, `operator()`, `equilibrium()`, `hamiltonian()`, `assume()`, `relaxation()`, `state()`, `evolution()`, `subplot()`.
