# examples/kinetics/relayed_hyperpol.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/examples/kinetics/relayed_hyperpol.m`
- Signature: `relayed_hyperpol()`
- Total lines: 101

## Purpose

Relayed NOE from hyperpolarized water to ALA-GLY dipeptide, generating Figure S7 from Christopher Pötzl

## Physical / mathematical content

- Chemical-kinetics examples. The files couple spin dynamics to exchange, pumping, or nonlinear reaction networks represented by kinetic generators in Liouville space.
- The relaxation model is Redfield-type perturbation theory: fluctuating interactions enter through correlation functions or spectral densities and generate a linear relaxation superoperator.

## Numerical / algorithmic content

- Time propagation is explicit. In Spinach this usually means repeated application of matrix exponentials or propagator factorizations to density operators or state vectors in Hilbert/Liouville/Fokker-Planck space.

## Code-derived implementation details

### Comment-guided execution stages

- Lines 10-11: Simulation timing parameters; implemented by `dt=0.125; nsteps=128`.
- Lines 13-14: Magnet field; implemented by `sys.magnet=16.4`.
- Lines 16-17: 30 protons in the system; implemented by `sys.isotopes=repelem({'1H'},30)`.
- Lines 19-20: Cartesian coordinates of pertinent protons; implemented by `inter.coordinates={[6.67 4.45 4.03]`.
- Lines 31-33: 20 water protons exist, but have no coordinates to prevent direct cross-relaxation from happening; implemented by `inter.coordinates=[inter.coordinates; repelem({[]},20)']`.
- Lines 35-37: Chemical shifts, all water at 4.5 ppm; implemented by `inter.zeeman.scalar={8.45 8.45 8.45 8.11 3.73 0.99 0.99 0.99 3.99 3.32}`.
- Lines 40-41: Relaxation theories; implemented by `inter.relaxation={'redfield','t1_t2'}`.
- Lines 47-48: Empirical relaxation at 0.1 Hz for water; implemented by `inter.r1_rates=num2cell([zeros(1,10) 0.1*ones(1,20)])`.
- Lines 51-53: Basis set, single-spin for water, up to three-spin orders for the molecule; implemented by `bas.formalism='sphten-liouv'`.
- Lines 59-60: Exchange flux matrix; implemented by `inter.chem.flux_rate=zeros(30,30)`.
- Lines 65-66: Spinach housekeeping; implemented by `spin_system=create(sys,inter)`.
- Lines 70-71: Build the Liouvillian; implemented by `H=hamiltonian(spin_system)`.
- Lines 76-77: Isotropic thermal equilibrium; implemented by `rho=equilibrium(spin_system)`.
- Lines 79-80: Polarise the water 100%; implemented by `Wz=state(spin_system,'Lz',11:20)`.
- Lines 83-84: Get detection states; implemented by `H_aliph=[6 7 8]; H_alpha=5`.
- Lines 88-90: Time evolution simulation; implemented by `result=evolution(spin_system,L,[HZ_aliph HZ_alpha], rho,dt,nsteps,'multichannel')`.
- Lines 92-93: Plotting; implemented by `time_axis=linspace(0,nsteps*dt,nsteps+1)`.

### Key state/data transformations

- Lines 11: computes `dt` using `dt=0.125; nsteps=128`.
- Lines 14: computes `sys.magnet` using `sys.magnet=16.4`.
- Lines 17: computes `sys.isotopes` using `sys.isotopes=repelem({'1H'},30)`.
- Lines 20: computes `inter.coordinates` using `inter.coordinates={[6.67 4.45 4.03]`.
- Lines 36-37: computes `inter.zeeman.scalar` using `inter.zeeman.scalar={8.45 8.45 8.45 8.11 3.73 0.99 0.99 0.99 3.99 3.32}`.
- Lines 41: computes `inter.relaxation` using `inter.relaxation={'redfield','t1_t2'}`.
- Lines 42: computes `inter.equilibrium` using `inter.equilibrium='dibari'`.
- Lines 43: computes `inter.rlx_keep` using `inter.rlx_keep='secular'`.
- Lines 44: computes `inter.tau_c` using `inter.tau_c={1.2e-10}`.
- Lines 45: computes `inter.temperature` using `inter.temperature=298`.
- Lines 48: computes `inter.r1_rates` using `inter.r1_rates=num2cell([zeros(1,10) 0.1*ones(1,20)])`.
- Lines 49: computes `inter.r2_rates` using `inter.r2_rates=num2cell([zeros(1,10) 0.1*ones(1,20)])`.
- Lines 53: computes `bas.formalism` using `bas.formalism='sphten-liouv'`.
- Lines 54: computes `bas.approximation` using `bas.approximation='IK-1'`.
- Lines 55: computes `bas.connectivity` using `bas.connectivity='full_tensors'`.
- Lines 56: computes `bas.space_level` using `bas.space_level=3`.
- Lines 57: computes `bas.level` using `bas.level=1`.
- Lines 60: computes `inter.chem.flux_rate` using `inter.chem.flux_rate=zeros(30,30)`.

## Implementation structure

- Relayed NOE from hyperpolarized water to ALA-GLY dipeptide,
- generating Figure S7 from
- Christopher Pötzl
- Simulation timing parameters
- Magnet field
- 30 protons in the system
- Cartesian coordinates of pertinent protons
- 20 water protons exist, but have no coordinates to
- prevent direct cross-relaxation from happening
- Chemical shifts, all water at 4.5 ppm
- Relaxation theories
- Empirical relaxation at 0.1 Hz for water

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `repelem()`, `num2cell()`, `create()`, `basis()`, `assume()`, `hamiltonian()`, `relaxation()`, `kinetics()`, `equilibrium()`, `state()`, `evolution()`, `figure()`, `scale_figure()`, `kxlabel()`, `kylabel()`, `klegend()`.
