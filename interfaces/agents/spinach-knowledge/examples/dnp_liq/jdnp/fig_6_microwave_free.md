# examples/dnp_liq/jdnp/fig_6_microwave_free.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/examples/dnp_liq/jdnp/fig_6_microwave_free.m`
- Signature: `fig_6_microwave_free()`
- Total lines: 134

## Purpose

A demonstration of Maria Grazia Concilio's microwave-free JDNP effect where a field ramp in combination with unequal relaxati- on rates of singlet-alpha and singlet-beta product states crea- tes nuclear magnetisation enhancement beyond the Boltzmann le- vel at both the starting and the final field. Calculation time: minutes

## Physical / mathematical content

- Liquid-state DNP examples. The main ingredients are electron-nuclear cross-relaxation, scalar or dipolar contact mechanisms, motional spectral densities, and field/frequency dependence of polarisation transfer.
- The relevant state manifold is the singlet/triplet decomposition, where permutation symmetry controls selection rules, relaxation susceptibility, and convertibility to ordinary magnetisation.

## Numerical / algorithmic content

- Time propagation is explicit. In Spinach this usually means repeated application of matrix exponentials or propagator factorizations to density operators or state vectors in Hilbert/Liouville/Fokker-Planck space.

## Code-derived implementation details

### Comment-guided execution stages

- Lines 14-15: Load the spin system; implemented by `[sys,inter,bas]=system_specification()`.
- Lines 17-18: Magnet fields; implemented by `start_field=14.09`.
- Lines 22-23: Match exchange coupling to the midpoint field; implemented by `inter.coupling.scalar{2,3}=match_field*(spin('E')+spin('1H'))/(2*pi)`.
- Lines 25-26: Increase viscosity; implemented by `inter.tau_c={2.2e-9}`.
- Lines 28-29: Get thermal equilibrium at starting field; implemented by `sys.magnet=start_field`.
- Lines 36-37: Set up a field ramp and time step; implemented by `field_grid=linspace(start_field,final_field,211); dt=1e-4`.
- Lines 39-40: Preallocate the trajectory; implemented by `traj=zeros(numel(rho_eq),numel(field_grid)+1,'like',1i)`.
- Lines 42-43: Set starting point; implemented by `traj(:,1)=rho_eq`.
- Lines 45-46: Run through the field ramp; implemented by `for n=1:numel(field_grid)`.
- Lines 48-49: Set the magneti field; implemented by `sys.magnet=field_grid(n)`.
- Lines 51-52: Run the housekeeping; implemented by `spin_system=create(sys,inter)`.
- Lines 55-56: Get the generators; implemented by `H=hamiltonian(assume(spin_system,'labframe'))`.
- Lines 59-61: Run the time evolution; implemented by `traj(:,n+1)=evolution(spin_system,H+1i*R,[],traj(:,n), dt,1,'final')`.
- Lines 65-66: Build component operators; implemented by `unit=unit_state(spin_system)`.
- Lines 79-80: Build alpha singlet and triplet states; implemented by `Sa=(unit/8)+(Nz/4)-(E1mE2p/4)-(E1pE2m/4)-(LzE2z/2)-(NzE1mE2p/2)-(NzE1pE2m/2)-NzE1zE2z`.
- Lines 85-86: Build beta singlet and triplet states; implemented by `Sb=(unit/8)-(Nz/4)-(E1mE2p/4)-(E1pE2m/4)-(LzE2z/2)+(NzE1mE2p/2)+(NzE1pE2m/2)+NzE1zE2z`.
- Lines 91-92: Build Nz-singlet and Nz-triplet states; implemented by `SNz=(Nz/4)-(NzE1mE2p/2)-(NzE1pE2m/2)-NzE1zE2z`.
- Lines 97-98: Detection states; implemented by `coils=[Tpa,Tpb,Tma,Tmb,T0a,T0b,Sa,Sb,SNz,TpNz,T0Nz,TmNz,E1z,E2z,Nz]`.

### Control flow inferred from the code

- Line 46: `for` loop over `n=1:numel(field_grid)`.

### Key state/data transformations

- Lines 15: computes `[sys,inter,bas]` using `[sys,inter,bas]=system_specification()`.
- Lines 18: computes `start_field` using `start_field=14.09`.
- Lines 19: computes `match_field` using `match_field=11.74`.
- Lines 20: computes `final_field` using `final_field=9.39`.
- Lines 23: computes `inter.coupling.scalar{2,3}` using `inter.coupling.scalar{2,3}=match_field*(spin('E')+spin('1H'))/(2*pi)`.
- Lines 26: computes `inter.tau_c` using `inter.tau_c={2.2e-9}`.
- Lines 29: computes `sys.magnet` using `sys.magnet=start_field`.
- Lines 30: computes `spin_system` using `spin_system=create(sys,inter)`.
- Lines 33: computes `H0` using `H0=hamiltonian(spin_system,'left')`.
- Lines 34: computes `rho_eq` using `rho_eq=equilibrium(spin_system,H0)`.
- Lines 37: computes `field_grid` using `field_grid=linspace(start_field,final_field,211); dt=1e-4`.
- Lines 40: computes `traj` using `traj=zeros(numel(rho_eq),numel(field_grid)+1,'like',1i)`.
- Lines 43: computes `traj(:,1)` using `traj(:,1)=rho_eq`.
- Lines 56: computes `H` using `H=hamiltonian(assume(spin_system,'labframe'))`.
- Lines 57: computes `R` using `R=relaxation(spin_system)`.
- Lines 60-61: computes `traj(:,n+1)` using `traj(:,n+1)=evolution(spin_system,H+1i*R,[],traj(:,n), dt,1,'final')`.
- Lines 66: computes `unit` using `unit=unit_state(spin_system)`.
- Lines 67: computes `Nz` using `Nz=state(spin_system,{'Lz'},{1})`.

## Implementation structure

- A demonstration of Maria Grazia Concilio's microwave-free JDNP
- effect where a field ramp in combination with unequal relaxati-
- on rates of singlet-alpha and singlet-beta product states crea-
- tes nuclear magnetisation enhancement beyond the Boltzmann le-
- vel at both the starting and the final field.
- Calculation time: minutes
- Load the spin system
- Magnet fields
- Match exchange coupling to the midpoint field
- Increase viscosity
- Get thermal equilibrium at starting field
- Set up a field ramp and time step

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `system_specification()`, `spin()`, `create()`, `basis()`, `assume()`, `hamiltonian()`, `equilibrium()`, `traj()`, `field_grid()`, `relaxation()`, `evolution()`, `unit_state()`, `state()`, `kfigure()`, `subplot()`, `scale_figure()`.
