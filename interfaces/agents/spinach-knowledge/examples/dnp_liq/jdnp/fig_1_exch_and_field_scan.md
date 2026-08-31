# examples/dnp_liq/jdnp/fig_1_exch_and_field_scan.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/examples/dnp_liq/jdnp/fig_1_exch_and_field_scan.m`
- Signature: `fig_1_exch_and_field_scan()`
- Total lines: 96

## Purpose

Matching condition plot for JDNP -proton polarisation at a particular time as a function of the external field and the inter-electron exchange coupling. Further details in: Calculation time: hours, line-by-line plotting

## Physical / mathematical content

- Liquid-state DNP examples. The main ingredients are electron-nuclear cross-relaxation, scalar or dipolar contact mechanisms, motional spectral densities, and field/frequency dependence of polarisation transfer.

## Numerical / algorithmic content

- Time propagation is explicit. In Spinach this usually means repeated application of matrix exponentials or propagator factorizations to density operators or state vectors in Hilbert/Liouville/Fokker-Planck space.
- The implementation explicitly addresses performance engineering through parallel or GPU execution, which matters because Spinach operators can become extremely large after basis expansion or powder/spatial lifting.

## Code-derived implementation details

### Comment-guided execution stages

- Lines 14-15: Load the spin system; implemented by `[sys,inter,bas,parameters]=system_specification()`.
- Lines 17-18: Experiment parameters; implemented by `parameters.mw_pwr=(2*pi*1e6)/2`.
- Lines 21-22: Field and coupling grids; implemented by `field_grid=linspace(0.25,3.0,256)`.
- Lines 25-26: Preallocate output array; implemented by `dnp=zeros([numel(field_grid) numel(exch_grid)])`.
- Lines 28-29: Create and scale the figure; implemented by `current_fig=kfigure(); scale_figure([1.0 0.65])`.
- Lines 31-32: Loop over the fields; implemented by `for n=1:numel(field_grid)`.
- Lines 34-35: Set the magnet field; implemented by `sys.magnet=field_grid(n)`.
- Lines 37-38: Trityl and free electron frequencies; implemented by `f_free=2*pi*g2freq(parameters.g_ref,sys.magnet)`.
- Lines 41-42: Microwave offset; implemented by `parameters.mw_off=f_trityl-f_free`.
- Lines 44-45: Parallel loop over exchange couplings; implemented by `parfor k=1:numel(exch_grid)`.
- Lines 47-48: Localisation; implemented by `locinter=inter`.
- Lines 50-51: Set the exchange coupling; implemented by `locinter.coupling.scalar=cell(3,3)`.
- Lines 54-55: Spinach housekeeping; implemented by `spin_system=create(sys,locinter)`.
- Lines 58-59: Relevant operators and states; implemented by `E_mw=operator(spin_system,'Lx','E')/2`.
- Lines 62-63: Isotropic thermal equilibrium; implemented by `rho_eq=equilibrium(spin_system)`.
- Lines 65-66: Detection state; implemented by `Nz=state(spin_system,'Lz','1H')`.
- Lines 68-69: Hamiltonian and relaxation superoperator; implemented by `H=hamiltonian(assume(spin_system,'esr'))`.
- Lines 72-73: Add microwave operator; implemented by `H=H+parameters.mw_pwr*E_mw`.

### Control flow inferred from the code

- Line 32: `for` loop over `n=1:numel(field_grid)`.
- Line 45: `parfor` loop over `k=1:numel(exch_grid)`.

### Key state/data transformations

- Lines 15: computes `[sys,inter,bas,parameters]` using `[sys,inter,bas,parameters]=system_specification()`.
- Lines 18: computes `parameters.mw_pwr` using `parameters.mw_pwr=(2*pi*1e6)/2`.
- Lines 19: computes `parameters.mw_dur` using `parameters.mw_dur=20e-3`.
- Lines 22: computes `field_grid` using `field_grid=linspace(0.25,3.0,256)`.
- Lines 23: computes `exch_grid` using `exch_grid=linspace(-100e9,100e9,256)`.
- Lines 26: computes `dnp` using `dnp=zeros([numel(field_grid) numel(exch_grid)])`.
- Lines 29: computes `current_fig` using `current_fig=kfigure(); scale_figure([1.0 0.65])`.
- Lines 35: computes `sys.magnet` using `sys.magnet=field_grid(n)`.
- Lines 38: computes `f_free` using `f_free=2*pi*g2freq(parameters.g_ref,sys.magnet)`.
- Lines 39: computes `f_trityl` using `f_trityl=2*pi*g2freq(parameters.g_trityl,sys.magnet)`.
- Lines 42: computes `parameters.mw_off` using `parameters.mw_off=f_trityl-f_free`.
- Lines 48: computes `locinter` using `locinter=inter`.
- Lines 51: computes `locinter.coupling.scalar` using `locinter.coupling.scalar=cell(3,3)`.
- Lines 52: computes `locinter.coupling.scalar{2,3}` using `locinter.coupling.scalar{2,3}=exch_grid(k)`.
- Lines 55: computes `spin_system` using `spin_system=create(sys,locinter)`.
- Lines 59: computes `E_mw` using `E_mw=operator(spin_system,'Lx','E')/2`.
- Lines 60: computes `Ez` using `Ez=operator(spin_system,'Lz','E')`.
- Lines 63: computes `rho_eq` using `rho_eq=equilibrium(spin_system)`.

## Implementation structure

- Matching condition plot for JDNP -proton polarisation at
- a particular time as a function of the external field and
- the inter-electron exchange coupling. Further details in:
- Calculation time: hours, line-by-line plotting
- Load the spin system
- Experiment parameters
- Field and coupling grids
- Preallocate output array
- Create and scale the figure
- Loop over the fields
- Set the magnet field
- Trityl and free electron frequencies

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `system_specification()`, `kfigure()`, `scale_figure()`, `field_grid()`, `g2freq()`, `exch_grid()`, `create()`, `basis()`, `operator()`, `equilibrium()`, `state()`, `hamiltonian()`, `assume()`, `relaxation()`, `evolution()`, `dnp()`.
