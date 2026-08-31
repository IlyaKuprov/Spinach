# tests/kernel/test_dynamic_pulse_utilities_deep.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/tests/kernel/test_dynamic_pulse_utilities_deep.m`
- Signature: `result=test_dynamic_pulse_utilities_deep()`
- Total lines: 219

## Purpose

Tests dynamic pulse utility paths. Syntax: result=test_dynamic_pulse_utilities_deep()

## Physical / mathematical content

- This file belongs to the `tests` part of Spinach. Its role should be read together with nearby files in the same directory, which usually share a common physical regime or infrastructure purpose.
- Orientation or trajectory averaging is performed numerically, so grid design, weights, and integration error control matter directly to accuracy and runtime.

## Numerical / algorithmic content

- Time propagation is explicit. In Spinach this usually means repeated application of matrix exponentials or propagator factorizations to density operators or state vectors in Hilbert/Liouville/Fokker-Planck space.
- Numerical integration over angles or geometry is part of the implementation, so point placement and weights are as important as the local Hamiltonian calculations.
- The file also defines local helper function(s): `local_liouv_system()`, `local_hilb_system()`, `grad_pulse_ref()`, `grad_sandw_ref()`, `local_delete()`. This usually means the public entry point is supported by tightly coupled validation or helper logic kept private to the file.

## Code-derived implementation details

### Comment-guided execution stages

- Lines 17-18: Announce the test target; implemented by `fprintf('TESTING: Dynamic pulse utility paths\n')`.
- Lines 20-23: State the pulse-utility target of the test; implemented by `result=new_test_result('kernel/dynamic_pulse_utilities_deep', 'Dynamic pulse utility paths', 'pulse utilities must produce finite dynamic outputs and match direct small-…`.
- Lines 25-26: Build a one-proton Liouville-space spin system with a carrier frequency; implemented by `spin_system=local_liouv_system(1.0)`.
- Lines 30-31: Compare grad_pulse with a direct small-matrix exponential reference; implemented by `rho_grad=grad_pulse(spin_system,L,rho,4.0,0.12,1.5e-4,0.8)`.
- Lines 38-39: Compare grad_sandw with a direct small-matrix exponential reference; implemented by `pulse_prop=propagator(spin_system,2*pi*50*operator(spin_system,'Ly',1),2.0e-4)`.
- Lines 47-48: Heterodyne a clean wall-clock carrier into the rotating frame; implemented by `dt=1e-4`.
- Lines 61-62: Exercise piecewise-constant and piecewise-linear RLC response models; implemented by `omega=2*pi*800`.
- Lines 79-80: Write a temporary Bruker pulse file and inspect stable JCAMP fields; implemented by `file_name=[tempname(fileparts(mfilename('fullpath'))) '.txt']`.
- Lines 93-94: Check finite-RF R-sequence compiler branches against direct exponentials; implemented by `spin_system_h=local_hilb_system()`.
- Lines 115-116: Check waveform basis variants on a compact odd grid; implemented by `basis_types={'sine_waves','cosine_waves','legendre'}`.
- Lines 119-120: Build and check each supported basis family; implemented by `basis_waves=wave_basis(basis_types{n},3,17)`.
- Lines 127-129: Check all supported pulse-shape variants against explicit formulae; implemented by `result=test_close(result,'pulse_shape rectangular deep',pulse_shape('rectangular',5),ones(1,5),1e-15,1e-15, 'rectangular pulse-shape samples are all unity')`.

### Control flow inferred from the code

- Line 117: `for` loop over `n=1:numel(basis_types)`.

### Key state/data transformations

- Lines 21-23: computes `result` using `result=new_test_result('kernel/dynamic_pulse_utilities_deep', 'Dynamic pulse utility paths', 'pulse utilities must produce finite dynamic outputs and match direct small-…`.
- Lines 26: computes `spin_system` using `spin_system=local_liouv_system(1.0)`.
- Lines 27: computes `L` using `L=0*operator(spin_system,'Lz',1)`.
- Lines 28: computes `rho` using `rho=state(spin_system,'Lx',1)`.
- Lines 31: computes `rho_grad` using `rho_grad=grad_pulse(spin_system,L,rho,4.0,0.12,1.5e-4,0.8)`.
- Lines 32: computes `rho_ref` using `rho_ref=grad_pulse_ref(spin_system,L,rho,4.0,0.12,1.5e-4,0.8)`.
- Lines 39: computes `pulse_prop` using `pulse_prop=propagator(spin_system,2*pi*50*operator(spin_system,'Ly',1),2.0e-4)`.
- Lines 40: computes `rho_sandw` using `rho_sandw=grad_sandw(spin_system,L,rho,pulse_prop,[3.0 -2.0],0.10,[1.0e-4 1.4e-4],[0.7 0.9])`.
- Lines 48: computes `dt` using `dt=1e-4`.
- Lines 49: computes `freq` using `freq=1000`.
- Lines 50: computes `time_grid` using `time_grid=dt*((1:4096)'-1)`.
- Lines 51: computes `signal` using `signal=cos(2*pi*freq*time_grid)`.
- Lines 52: computes `[X_het,Y_het]` using `[X_het,Y_het]=heterodyne(dt,signal,freq)`.
- Lines 53: computes `steady` using `steady=200:(numel(signal)-200)`.
- Lines 62: computes `omega` using `omega=2*pi*800`.
- Lines 63: computes `Q` using `Q=8`.
- Lines 64: computes `dt_user` using `dt_user=1e-3`.
- Lines 65: computes `[X_pwc,Y_pwc,dt_pwc]` using `[X_pwc,Y_pwc,dt_pwc]=restrans([1;2;1;0],[0;0.5;0;-0.5],dt_user,omega,Q,'pwc',2)`.

### Local helper functions

- Line 142: `local_liouv_system()` — `function spin_system=local_liouv_system(magnet)`. Specify the spin system
  - Representative operation: `sys.magnet=magnet`.
  - Representative operation: `sys.isotopes={'1H'}`.
- Line 159: `local_hilb_system()` — `function spin_system=local_hilb_system()`. Specify the spin system
  - Representative operation: `sys.magnet=0`.
  - Representative operation: `sys.isotopes={'1H'}`.
- Line 176: `grad_pulse_ref()` — `function rho=grad_pulse_ref(spin_system,L,rho,g_amp,s_len,g_dur,s_fac)`. Build the effective gradient operator
  - Representative operation: `G=1e-4*s_fac*g_amp*s_len*g_dur*carrier(spin_system,'all')/spin_system.inter.magnet`.
  - Representative operation: `rho=expm(-1i*full(L)*g_dur)*rho`.
- Line 191: `grad_sandw_ref()` — `function rho=grad_sandw_ref(spin_system,L,rho,P,g_amps,s_len,g_durs,s_facs)`. Build the effective gradient operators
  - Representative operation: `R=carrier(spin_system,'all')/spin_system.inter.magnet`.
  - Representative operation: `G1=1e-4*s_facs(1)*g_amps(1)*s_len*g_durs(1)*R`.
- Line 210: `local_delete()` — `function local_delete(file_name)`. Delete the file if it exists
  - Representative operation: `if exist(file_name,'file')==2`.
  - Representative operation: `delete(file_name)`.

## Outputs

- result -regression test result with explanatory messages
- The test checks gradient-pulse dynamics, heterodyne filtering, RLC
- response transforms, Bruker pulse-file writing, finite-RF R-sequence
- compilation, waveform basis variants, and pulse-shape variants.

## Implementation structure

- Tests dynamic pulse utility paths. Syntax:
- result=test_dynamic_pulse_utilities_deep()
- result -regression test result with explanatory messages
- The test checks gradient-pulse dynamics, heterodyne filtering, RLC
- response transforms, Bruker pulse-file writing, finite-RF R-sequence
- compilation, waveform basis variants, and pulse-shape variants.
- Announce the test target
- State the pulse-utility target of the test
- Build a one-proton Liouville-space spin system with a carrier frequency
- Compare grad_pulse with a direct small-matrix exponential reference
- Compare grad_sandw with a direct small-matrix exponential reference
- Heterodyne a clean wall-clock carrier into the rotating frame

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `new_test_result()`, `local_liouv_system()`, `operator()`, `state()`, `grad_pulse()`, `grad_pulse_ref()`, `test_close()`, `test_true()`, `propagator()`, `grad_sandw()`, `grad_sandw_ref()`, `heterodyne()`, `X_het()`, `Y_het()`, `all()`, `restrans()`.
