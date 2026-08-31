# examples/dnp_liq/jdnp/fig_3_time_dep_top_row.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/examples/dnp_liq/jdnp/fig_3_time_dep_top_row.m`
- Signature: `fig_3_time_dep_top_row()`
- Total lines: 87

## Purpose

A demonstration that the JDNP effect vanishes when the second electron is removed from the system. Proton polarisation as a function of time for specific external fields is plotted. See also the "bot row" simulation where both electrons are active and the JDNP enhancement is present. Further details in: Calculation time: seconds, line-by-line plotting

## Physical / mathematical content

- Liquid-state DNP examples. The main ingredients are electron-nuclear cross-relaxation, scalar or dipolar contact mechanisms, motional spectral densities, and field/frequency dependence of polarisation transfer.
- The relaxation model is Redfield-type perturbation theory: fluctuating interactions enter through correlation functions or spectral densities and generate a linear relaxation superoperator.

## Numerical / algorithmic content

- Time propagation is explicit. In Spinach this usually means repeated application of matrix exponentials or propagator factorizations to density operators or state vectors in Hilbert/Liouville/Fokker-Planck space.

## Code-derived implementation details

### Comment-guided execution stages

- Lines 16-17: Load the spin system; implemented by `[sys,inter,bas,parameters]=system_specification()`.
- Lines 19-20: Kill the second electron; implemented by `sys.isotopes=sys.isotopes([1 2])`.
- Lines 27-28: Experiment parameters; implemented by `parameters.mw_pwr=2*pi*250e3`.
- Lines 32-33: Magnetic field grid, Tesla; implemented by `field_grid=[0.034 0.34 3.4]`.
- Lines 35-36: Get a figure going; implemented by `kfigure(); scale_figure([1.75 0.5])`.
- Lines 38-39: Loop over the fields; implemented by `for n=1:numel(field_grid)`.
- Lines 41-42: Set magnet field; implemented by `sys.magnet=field_grid(n)`.
- Lines 44-45: Set microwave offset frequency; implemented by `f_free=g2freq(parameters.g_ref,sys.magnet)`.
- Lines 49-50: Spinach housekeeping; implemented by `spin_system=create(sys,inter)`.
- Lines 53-54: Electron operators; implemented by `Ex=operator(spin_system,'Lx','E')`.
- Lines 57-58: Isotropic thermal equilibrium; implemented by `rho_eq=equilibrium(spin_system)`.
- Lines 60-61: Hamiltonian and relaxation superoperator; implemented by `H=hamiltonian(assume(spin_system,'esr'))`.
- Lines 64-65: Add microwave terms to the Hamiltonian; implemented by `H=H+parameters.mw_pwr*Ex+parameters.mw_off*Ez`.
- Lines 67-68: Detection state: Lz on the proton; implemented by `coil=state(spin_system,'Lz','1H')`.
- Lines 70-72: Run the time evolution and normalise to thermal equilibrium; implemented by `answer=evolution(spin_system,H+1i*R,coil,rho_eq,parameters.t_step, parameters.nsteps,'observable')`.
- Lines 75-77: Time axis generation; implemented by `t_axis=linspace(0,parameters.t_step*parameters.nsteps, parameters.nsteps+1)`.
- Lines 79-80: Plotting; implemented by `subplot(1,3,n); plot(t_axis,answer); kgrid; box on`.

### Control flow inferred from the code

- Line 39: `for` loop over `n=1:numel(field_grid)`.

### Key state/data transformations

- Lines 17: computes `[sys,inter,bas,parameters]` using `[sys,inter,bas,parameters]=system_specification()`.
- Lines 20: computes `sys.isotopes` using `sys.isotopes=sys.isotopes([1 2])`.
- Lines 21: computes `inter.zeeman.matrix` using `inter.zeeman.matrix=inter.zeeman.matrix([1 2])`.
- Lines 22: computes `inter.coordinates` using `inter.coordinates=inter.coordinates([1 2])`.
- Lines 23: computes `inter.coupling.scalar` using `inter.coupling.scalar=cell(2,2)`.
- Lines 24: computes `inter` using `inter=rmfield(inter,{'srfk_tau_c','srfk_mdepth'})`.
- Lines 25: computes `inter.relaxation` using `inter.relaxation={'redfield'}`.
- Lines 28: computes `parameters.mw_pwr` using `parameters.mw_pwr=2*pi*250e3`.
- Lines 29: computes `parameters.t_step` using `parameters.t_step=1e-3`.
- Lines 30: computes `parameters.nsteps` using `parameters.nsteps=300`.
- Lines 33: computes `field_grid` using `field_grid=[0.034 0.34 3.4]`.
- Lines 42: computes `sys.magnet` using `sys.magnet=field_grid(n)`.
- Lines 45: computes `f_free` using `f_free=g2freq(parameters.g_ref,sys.magnet)`.
- Lines 46: computes `f_trityl` using `f_trityl=g2freq(parameters.g_trityl,sys.magnet)`.
- Lines 47: computes `parameters.mw_off` using `parameters.mw_off=2*pi*(f_trityl-f_free)`.
- Lines 50: computes `spin_system` using `spin_system=create(sys,inter)`.
- Lines 54: computes `Ex` using `Ex=operator(spin_system,'Lx','E')`.
- Lines 55: computes `Ez` using `Ez=operator(spin_system,'Lz','E')`.

## Implementation structure

- A demonstration that the JDNP effect vanishes when the second
- electron is removed from the system. Proton polarisation as a
- function of time for specific external fields is plotted. See
- also the "bot row" simulation where both electrons are active
- and the JDNP enhancement is present. Further details in:
- Calculation time: seconds, line-by-line plotting
- Load the spin system
- Kill the second electron
- Experiment parameters
- Magnetic field grid, Tesla
- Get a figure going
- Loop over the fields

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `system_specification()`, `rmfield()`, `kfigure()`, `scale_figure()`, `field_grid()`, `g2freq()`, `create()`, `basis()`, `operator()`, `equilibrium()`, `hamiltonian()`, `assume()`, `relaxation()`, `state()`, `evolution()`, `subplot()`.
