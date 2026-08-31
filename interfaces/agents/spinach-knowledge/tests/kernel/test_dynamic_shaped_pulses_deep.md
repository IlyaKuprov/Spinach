# tests/kernel/test_dynamic_shaped_pulses_deep.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/tests/kernel/test_dynamic_shaped_pulses_deep.m`
- Signature: `result=test_dynamic_shaped_pulses_deep()`
- Total lines: 140

## Purpose

Tests dynamic shaped-pulse propagation paths. Syntax: result=test_dynamic_shaped_pulses_deep()

## Physical / mathematical content

- This file belongs to the `tests` part of Spinach. Its role should be read together with nearby files in the same directory, which usually share a common physical regime or infrastructure purpose.
- The file uses a Fokker-Planck-style enlarged state space in which spatial or orientational coordinates are promoted to extra dimensions and coupled to spin dynamics through differential operators.
- Propagation is accelerated with a Krylov-subspace method, replacing direct matrix exponentiation by projection into a much smaller Arnoldi/Lanczos-type subspace.
- Orientation or trajectory averaging is performed numerically, so grid design, weights, and integration error control matter directly to accuracy and runtime.

## Numerical / algorithmic content

- Time propagation is explicit. In Spinach this usually means repeated application of matrix exponentials or propagator factorizations to density operators or state vectors in Hilbert/Liouville/Fokker-Planck space.
- A Krylov-subspace or Arnoldi construction is used to avoid forming or exponentiating very large dense propagators directly.
- Numerical integration over angles or geometry is part of the implementation, so point placement and weights are as important as the local Hamiltonian calculations.
- The file also defines local helper function(s): `local_liouv_system()`. This usually means the public entry point is supported by tightly coupled validation or helper logic kept private to the file.

## Code-derived implementation details

### Comment-guided execution stages

- Lines 17-18: Announce the test target; implemented by `fprintf('TESTING: Dynamic shaped pulse propagation paths\n')`.
- Lines 20-23: State the dynamic shaped-pulse target of the test; implemented by `result=new_test_result('kernel/dynamic_shaped_pulses_deep', 'Dynamic shaped pulse propagation paths', 'shaped-pulse helpers must reduce to exact constant-generator propa…`.
- Lines 25-26: Build a one-proton Liouville-space spin system; implemented by `spin_system=local_liouv_system(0)`.
- Lines 28-29: Get compact Liouville-space controls; implemented by `Lx=operator(spin_system,'Lx',1)`.
- Lines 34-35: Define a constant Cartesian RF generator over two slices; implemented by `slice_durs=[8e-5 13e-5]`.
- Lines 43-46: Exercise the Krylov piecewise-constant path; implemented by `[rho_obs,traj_obs]=shaped_pulse_xy(spin_system,drift,{Lx,Ly}, {[amp_x amp_x],[amp_y amp_y]}, slice_durs,rho,'expv-pwc')`.
- Lines 54-57: Exercise the Krylov piecewise-linear path; implemented by `[rho_obs,traj_obs]=shaped_pulse_xy(spin_system,drift,{Lx,Ly}, {amp_x*ones(1,3),amp_y*ones(1,3)}, slice_durs,rho,'expv-pwl')`.
- Lines 63-64: Exercise explicit exponential product quadratures and propagator output; implemented by `methods={'expm-pwc','expm-pwl','evol-pwc','evol-pwl'}`.
- Lines 67-68: Select the matching amplitude table shape; implemented by `if strcmp(methods{n}(6:8),'pwc')`.
- Lines 74-76: Run the shaped pulse with propagator output; implemented by `[rho_obs,traj_obs,prop_obs]=shaped_pulse_xy(spin_system,drift,{Lx,Ly}, amplitudes,slice_durs,rho,methods{n})`.
- Lines 78-80: Check final state and propagator consistency; implemented by `result=test_close(result,['shaped_pulse_xy ' methods{n} ' final'],rho_obs,ref_rho,1e-10,1e-10, 'matrix-exponential and evolution paths reduce to exact constant-generator…`.
- Lines 87-88: Define a constant amplitude-frequency pulse; implemented by `rf_phi=pi/7`.
- Lines 95-96: Exercise all shaped_pulse_af propagation methods; implemented by `methods={'expv','expm','evolution'}`.
- Lines 99-100: Run the Fokker-Planck shaped pulse; implemented by `if strcmp(methods{n},'expm')`.
- Lines 114-116: Check the folded final state and trajectory; implemented by `result=test_close(result,['shaped_pulse_af ' methods{n} ' final'],rho_obs,ref_rho,1e-10,1e-10, 'constant amplitude-frequency controls fold back to exact Cartesian propag…`.

### Control flow inferred from the code

- Line 65: `for` loop over `n=1:numel(methods)`.
- Line 68: conditional branch on `strcmp(methods{n}(6:8),'pwc')`.
- Line 97: `for` loop over `n=1:numel(methods)`.
- Line 100: conditional branch on `strcmp(methods{n},'expm')`.

### Key state/data transformations

- Lines 21-23: computes `result` using `result=new_test_result('kernel/dynamic_shaped_pulses_deep', 'Dynamic shaped pulse propagation paths', 'shaped-pulse helpers must reduce to exact constant-generator propa…`.
- Lines 26: computes `spin_system` using `spin_system=local_liouv_system(0)`.
- Lines 29: computes `Lx` using `Lx=operator(spin_system,'Lx',1)`.
- Lines 30: computes `Ly` using `Ly=operator(spin_system,'Ly',1)`.
- Lines 31: computes `Lz` using `Lz=operator(spin_system,'Lz',1)`.
- Lines 32: computes `rho` using `rho=state(spin_system,'Lz',1)`.
- Lines 35: computes `slice_durs` using `slice_durs=[8e-5 13e-5]`.
- Lines 36: computes `amp_x` using `amp_x=2*pi*430`.
- Lines 37: computes `amp_y` using `amp_y=-2*pi*170`.
- Lines 38: computes `drift` using `drift=2*pi*35*Lz`.
- Lines 39: computes `generator` using `generator=drift+amp_x*Lx+amp_y*Ly`.
- Lines 40: computes `ref_rho` using `ref_rho=step(spin_system,generator,rho,sum(slice_durs))`.
- Lines 41: computes `ref_prop` using `ref_prop=propagator(spin_system,generator,sum(slice_durs))`.
- Lines 44-46: computes `[rho_obs,traj_obs]` using `[rho_obs,traj_obs]=shaped_pulse_xy(spin_system,drift,{Lx,Ly}, {[amp_x amp_x],[amp_y amp_y]}, slice_durs,rho,'expv-pwc')`.
- Lines 64: computes `methods` using `methods={'expm-pwc','expm-pwl','evol-pwc','evol-pwl'}`.
- Lines 69: computes `amplitudes` using `amplitudes={[amp_x amp_x],[amp_y amp_y]}`.
- Lines 75-76: computes `[rho_obs,traj_obs,prop_obs]` using `[rho_obs,traj_obs,prop_obs]=shaped_pulse_xy(spin_system,drift,{Lx,Ly}, amplitudes,slice_durs,rho,methods{n})`.
- Lines 88: computes `rf_phi` using `rf_phi=pi/7`.

### Local helper functions

- Line 124: `local_liouv_system()` — `function spin_system=local_liouv_system(magnet)`. Specify the spin system
  - Representative operation: `sys.magnet=magnet`.
  - Representative operation: `sys.isotopes={'1H'}`.

## Outputs

- result -regression test result with explanatory messages
- The test checks constant-generator reductions of shaped_pulse_xy
- product-quadrature/method combinations and shaped_pulse_af Fokker-Planck
- propagation paths.

## Implementation structure

- Tests dynamic shaped-pulse propagation paths. Syntax:
- result=test_dynamic_shaped_pulses_deep()
- result -regression test result with explanatory messages
- The test checks constant-generator reductions of shaped_pulse_xy
- product-quadrature/method combinations and shaped_pulse_af Fokker-Planck
- propagation paths.
- Announce the test target
- State the dynamic shaped-pulse target of the test
- Build a one-proton Liouville-space spin system
- Get compact Liouville-space controls
- Define a constant Cartesian RF generator over two slices
- Exercise the Krylov piecewise-constant path

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `new_test_result()`, `local_liouv_system()`, `operator()`, `state()`, `step()`, `propagator()`, `shaped_pulse_xy()`, `test_close()`, `strcmp()`, `shaped_pulse_af()`, `test_spin_system()`.
