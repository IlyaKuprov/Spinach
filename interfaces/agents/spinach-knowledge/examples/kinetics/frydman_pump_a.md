# examples/kinetics/frydman_pump_a.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/examples/kinetics/frydman_pump_a.m`
- Signature: `frydman_pump_a()`
- Total lines: 132

## Purpose

Lucio Frydman's water exchange based spin-lock pump, Figure 2 from https://doi.org/10.1016/j.jmr.2021.107083 Calculation time: seconds.

## Physical / mathematical content

- Chemical-kinetics examples. The files couple spin dynamics to exchange, pumping, or nonlinear reaction networks represented by kinetic generators in Liouville space.
- Propagation is accelerated with a Krylov-subspace method, replacing direct matrix exponentiation by projection into a much smaller Arnoldi/Lanczos-type subspace.

## Numerical / algorithmic content

- Time propagation is explicit. In Spinach this usually means repeated application of matrix exponentials or propagator factorizations to density operators or state vectors in Hilbert/Liouville/Fokker-Planck space.
- A Krylov-subspace or Arnoldi construction is used to avoid forming or exponentiating very large dense propagators directly.
- The file also defines local helper function(s): `frydman_pump()`. This usually means the public entry point is supported by tightly coupled validation or helper logic kept private to the file.

## Code-derived implementation details

### Comment-guided execution stages

- Lines 13-14: Number of water protons; implemented by `n_water_protons=20`.
- Lines 16-17: Magnet field; implemented by `sys.magnet=11.7`.
- Lines 19-20: Core spin system; implemented by `sys.isotopes={'1H','15N','13C','13C'}`.
- Lines 22-23: Add water protons; implemented by `sys.isotopes(5:(5+n_water_protons-1))={'1H'}`.
- Lines 25-26: Chemical shifts; implemented by `inter.zeeman.scalar={0.0,0.0,0.0,0.0}`.
- Lines 29-30: Scalar couplings, Hz; implemented by `inter.coupling.scalar=cell(n_water_protons+4)`.
- Lines 33-34: Estimated relaxation times, seconds; implemented by `T1H=0.2722; T1N=0.8; T1Ca=2; T1CO=2; T1W=0.2994`.
- Lines 37-38: Relaxation theory; implemented by `inter.relaxation={'t1_t2'}`.
- Lines 45-46: Basis set; implemented by `bas.formalism='sphten-liouv'`.
- Lines 52-53: Exchange rates, Hz; implemented by `nh_wt_exch_rate=10`.
- Lines 56-57: Exchange rate matrix; implemented by `inter.chem.flux_rate=zeros(n_water_protons+4)`.
- Lines 64-65: Spinach housekeeping; implemented by `spin_system=create(sys,inter)`.
- Lines 68-69: Simulation parameters; implemented by `parameters.spins={'1H','15N','13C'}`.
- Lines 73-74: Get system trajectory (sequence is below); implemented by `traj=liquid(spin_system,@frydman_pump,parameters,'nmr')`.
- Lines 76-77: Observables: peptide bond N-H; implemented by `Hz=state(spin_system,{'Lz'},{1})`.
- Lines 82-83: Project out the observables; implemented by `Hz=real(Hz'*traj); Hx=real(Hx'*traj)`.
- Lines 86-87: Plotting; implemented by `kfigure(); scale_figure([1.0 1.5])`.

### Key state/data transformations

- Lines 14: computes `n_water_protons` using `n_water_protons=20`.
- Lines 17: computes `sys.magnet` using `sys.magnet=11.7`.
- Lines 20: computes `sys.isotopes` using `sys.isotopes={'1H','15N','13C','13C'}`.
- Lines 23: computes `sys.isotopes(5:(5+n_water_protons-1))` using `sys.isotopes(5:(5+n_water_protons-1))={'1H'}`.
- Lines 26: computes `inter.zeeman.scalar` using `inter.zeeman.scalar={0.0,0.0,0.0,0.0}`.
- Lines 27: computes `inter.zeeman.scalar(5:(5+n_water_protons-1))` using `inter.zeeman.scalar(5:(5+n_water_protons-1))={0.0}`.
- Lines 30: computes `inter.coupling.scalar` using `inter.coupling.scalar=cell(n_water_protons+4)`.
- Lines 31: computes `inter.coupling.scalar{1,2}` using `inter.coupling.scalar{1,2}=-45`.
- Lines 34: computes `T1H` using `T1H=0.2722; T1N=0.8; T1Ca=2; T1CO=2; T1W=0.2994`.
- Lines 35: computes `T2H` using `T2H=0.2722; T2N=0.8; T2Ca=2; T2CO=2; T2W=0.2994`.
- Lines 38: computes `inter.relaxation` using `inter.relaxation={'t1_t2'}`.
- Lines 39: computes `inter.r1_rates` using `inter.r1_rates=num2cell(1./[T1H T1N T1Ca T1CO T1W*ones(1,n_water_protons)])`.
- Lines 40: computes `inter.r2_rates` using `inter.r2_rates=num2cell(1./[T2H T2N T2Ca T2CO T2W*ones(1,n_water_protons)])`.
- Lines 41: computes `inter.rlx_keep` using `inter.rlx_keep='diagonal'`.
- Lines 42: computes `inter.equilibrium` using `inter.equilibrium='IME'`.
- Lines 43: computes `inter.temperature` using `inter.temperature=298`.
- Lines 46: computes `bas.formalism` using `bas.formalism='sphten-liouv'`.
- Lines 47: computes `bas.approximation` using `bas.approximation='IK-1'`.

### Local helper functions

- Line 98: `frydman_pump()` — `function traj=frydman_pump(spin_system,parameters,~,R,K)`. Get pulse operators
  - Representative operation: `Hy=operator(spin_system,'Ly',parameters.spins{1})`.
  - Representative operation: `Ny=operator(spin_system,'Ly',parameters.spins{2})`.

## Implementation structure

- Lucio Frydman's water exchange based spin-lock pump, Figure 2
- from https://doi.org/10.1016/j.jmr.2021.107083
- Calculation time: seconds.
- Number of water protons
- Magnet field
- Core spin system
- Add water protons
- Chemical shifts
- Scalar couplings, Hz
- Estimated relaxation times, seconds
- Relaxation theory
- Basis set

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `num2cell()`, `create()`, `basis()`, `liquid()`, `state()`, `kfigure()`, `scale_figure()`, `subplot()`, `klegend()`, `kylabel()`, `kxlabel()`, `frydman_pump()`, `operator()`, `equilibrium()`, `dictum()`, `hamiltonian()`.
