# examples/dnp_liq/jdnp/fig_4_state_amplitudes.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/examples/dnp_liq/jdnp/fig_4_state_amplitudes.m`
- Signature: `fig_4_state_amplitudes()`
- Total lines: 126

## Purpose

Time evolution of the individual states in the basis set built from the singlet-triplet basis on the two electrons and Carte- sian spin operator basis on the nucleus, demonstrating the im- balance between singlet-alpha and singlet-beta subspace. This leads to trainsient nuclear spin polarisation. Details in: Calculation time: seconds

## Physical / mathematical content

- Liquid-state DNP examples. The main ingredients are electron-nuclear cross-relaxation, scalar or dipolar contact mechanisms, motional spectral densities, and field/frequency dependence of polarisation transfer.
- The relevant state manifold is the singlet/triplet decomposition, where permutation symmetry controls selection rules, relaxation susceptibility, and convertibility to ordinary magnetisation.

## Numerical / algorithmic content

- Time propagation is explicit. In Spinach this usually means repeated application of matrix exponentials or propagator factorizations to density operators or state vectors in Hilbert/Liouville/Fokker-Planck space.

## Code-derived implementation details

### Comment-guided execution stages

- Lines 16-17: Load the spin system; implemented by `[sys,inter,bas,parameters]=system_specification()`.
- Lines 19-20: Experiment parameters; implemented by `parameters.mw_pwr=2*pi*250e3`.
- Lines 24-25: Set the magnet; implemented by `sys.magnet=14.08`.
- Lines 27-28: Set microwave offset frequency; implemented by `f_free=g2freq(parameters.g_ref,sys.magnet)`.
- Lines 32-33: Set the exchange coupling; implemented by `electron_zeeman_tensor=g2freq(inter.zeeman.matrix{2},sys.magnet)`.
- Lines 38-39: Spinach housekeeping; implemented by `spin_system=create(sys,inter)`.
- Lines 42-43: Electron operators; implemented by `Ex=operator(spin_system,'Lx','E')`.
- Lines 46-47: Isotropic thermal equilibrium; implemented by `rho_eq=equilibrium(spin_system)`.
- Lines 49-50: Hamiltonian and relaxation superoperator; implemented by `H=hamiltonian(assume(spin_system,'esr'))`.
- Lines 53-54: Add microwave terms to the Hamiltonian; implemented by `H=H+parameters.mw_pwr*Ex+parameters.mw_off*Ez`.
- Lines 56-57: Build component operators; implemented by `unit=unit_state(spin_system)`.
- Lines 70-71: Build alpha singlet and triplet states; implemented by `Sa=(unit/8)+(Nz/4)-(E1mE2p/4)-(E1pE2m/4)-(LzE2z/2)-(NzE1mE2p/2)-(NzE1pE2m/2)-NzE1zE2z`.
- Lines 76-77: Build beta singlet and triplet states; implemented by `Sb=(unit/8)-(Nz/4)-(E1mE2p/4)-(E1pE2m/4)-(LzE2z/2)+(NzE1mE2p/2)+(NzE1pE2m/2)+NzE1zE2z`.
- Lines 82-83: Build Nz-singlet and Nz-triplet states; implemented by `SNz=(Nz/4)-(NzE1mE2p/2)-(NzE1pE2m/2)-NzE1zE2z`.
- Lines 88-89: Detection states; implemented by `coils=[Tpa,Tpb,Tma,Tmb,T0a,T0b,Sa,Sb,SNz,TpNz,T0Nz,TmNz,E1z,E2z,Nz]`.
- Lines 91-93: Run the time evolution; implemented by `answer=evolution(spin_system,H+1i*R,coils,rho_eq,parameters.t_step, parameters.nsteps,'multichannel')`.
- Lines 95-96: Time axis generation; implemented by `t_axis=linspace(0,parameters.t_step*parameters.nsteps,parameters.nsteps+1)`.
- Lines 98-99: Plotting; implemented by `kfigure(); scale_figure([1.75 0.5]); subplot(1,3,1)`.

### Key state/data transformations

- Lines 17: computes `[sys,inter,bas,parameters]` using `[sys,inter,bas,parameters]=system_specification()`.
- Lines 20: computes `parameters.mw_pwr` using `parameters.mw_pwr=2*pi*250e3`.
- Lines 21: computes `parameters.t_step` using `parameters.t_step=1e-3`.
- Lines 22: computes `parameters.nsteps` using `parameters.nsteps=700`.
- Lines 25: computes `sys.magnet` using `sys.magnet=14.08`.
- Lines 28: computes `f_free` using `f_free=g2freq(parameters.g_ref,sys.magnet)`.
- Lines 29: computes `f_trityl` using `f_trityl=g2freq(parameters.g_trityl,sys.magnet)`.
- Lines 30: computes `parameters.mw_off` using `parameters.mw_off=2*pi*(f_trityl-f_free)`.
- Lines 33: computes `electron_zeeman_tensor` using `electron_zeeman_tensor=g2freq(inter.zeeman.matrix{2},sys.magnet)`.
- Lines 34: computes `electron_zeeman_iso` using `electron_zeeman_iso=mean(diag(electron_zeeman_tensor))`.
- Lines 35: computes `proton_zeeman_iso` using `proton_zeeman_iso=sys.magnet*spin('1H')/(2*pi)`.
- Lines 36: computes `inter.coupling.scalar{2,3}` using `inter.coupling.scalar{2,3}=electron_zeeman_iso+proton_zeeman_iso`.
- Lines 39: computes `spin_system` using `spin_system=create(sys,inter)`.
- Lines 43: computes `Ex` using `Ex=operator(spin_system,'Lx','E')`.
- Lines 44: computes `Ez` using `Ez=operator(spin_system,'Lz','E')`.
- Lines 47: computes `rho_eq` using `rho_eq=equilibrium(spin_system)`.
- Lines 50: computes `H` using `H=hamiltonian(assume(spin_system,'esr'))`.
- Lines 51: computes `R` using `R=relaxation(spin_system)`.

## Implementation structure

- Time evolution of the individual states in the basis set built
- from the singlet-triplet basis on the two electrons and Carte-
- sian spin operator basis on the nucleus, demonstrating the im-
- balance between singlet-alpha and singlet-beta subspace. This
- leads to trainsient nuclear spin polarisation. Details in:
- Calculation time: seconds
- Load the spin system
- Experiment parameters
- Set the magnet
- Set microwave offset frequency
- Set the exchange coupling
- Spinach housekeeping

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `system_specification()`, `g2freq()`, `spin()`, `create()`, `basis()`, `operator()`, `equilibrium()`, `hamiltonian()`, `assume()`, `relaxation()`, `unit_state()`, `state()`, `evolution()`, `kfigure()`, `scale_figure()`, `subplot()`.
