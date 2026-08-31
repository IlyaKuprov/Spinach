# experiments/rapidscan.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/experiments/rapidscan.m`
- Signature: `[b_axis,spectrum]=rapidscan(spin_system,parameters)`
- Total lines: 120

## Purpose

Time-domain rapid field scan ESR experiment, Eatons style. Syntax: [b_axis,spectrum]=rapidscan(spin_system,parameters)

## Physical / mathematical content

- This file belongs to the `experiments` part of Spinach. Its role should be read together with nearby files in the same directory, which usually share a common physical regime or infrastructure purpose.

## Numerical / algorithmic content

- Time propagation is explicit. In Spinach this usually means repeated application of matrix exponentials or propagator factorizations to density operators or state vectors in Hilbert/Liouville/Fokker-Planck space.
- The file contains an explicit `grumble(...)` validator, which is Spinach convention for front-loading dimension, type, and regime checks before expensive linear-algebra work begins.
- The file also defines local helper function(s): `grumble()`. This usually means the public entry point is supported by tightly coupled validation or helper logic kept private to the file.

## Code-derived implementation details

### Comment-guided execution stages

- Lines 34-35: Check consistency; implemented by `grumble(parameters)`.
- Lines 37-38: Get the microwave Hamiltonian; implemented by `Ep=operator(spin_system,'L+','E')`.
- Lines 41-42: Isotropic thermal equilibrium; implemented by `rho=equilibrium(spin_system)`.
- Lines 44-45: Relevant operators; implemented by `Hz=hamiltonian(assume(spin_system,'labframe','zeeman'))`.
- Lines 49-50: Put the drift Liouvillian together; implemented by `L=Hz+Hc; L=(L+L')/2`.
- Lines 52-53: Go into the microwave rotating frame; implemented by `spin_system=assume(spin_system,'labframe')`.
- Lines 57-58: Add microwave and relaxation terms; implemented by `L=L+H_mw+1i*R`.
- Lines 60-61: Get the detection state; implemented by `coil=state(spin_system,'L+','E')`.
- Lines 63-65: Compute the waveform and the axis; implemented by `waveform=linspace(parameters.sweep(1), parameters.sweep(2),parameters.nsteps)`.
- Lines 68-69: Normalize the Zeeman Hamiltonian; implemented by `Hz=Hz/spin_system.inter.magnet`.
- Lines 71-72: Preallocate the answer; implemented by `spectrum=zeros([parameters.nsteps 1],'like',1i)`.
- Lines 74-75: Propagate the system; implemented by `for n=1:parameters.nsteps`.

### Control flow inferred from the code

- Line 75: `for` loop over `n=1:parameters.nsteps`.

### Key state/data transformations

- Lines 38: computes `Ep` using `Ep=operator(spin_system,'L+','E')`.
- Lines 39: computes `H_mw` using `H_mw=-parameters.mw_pwr*(Ep-Ep')/2i`.
- Lines 42: computes `rho` using `rho=equilibrium(spin_system)`.
- Lines 45: computes `Hz` using `Hz=hamiltonian(assume(spin_system,'labframe','zeeman'))`.
- Lines 46: computes `Hc` using `Hc=hamiltonian(assume(spin_system,'labframe','couplings'))`.
- Lines 47: computes `R` using `R=relaxation(spin_system)`.
- Lines 50: computes `L` using `L=Hz+Hc; L=(L+L')/2`.
- Lines 53: computes `spin_system` using `spin_system=assume(spin_system,'labframe')`.
- Lines 54: computes `C` using `C=carrier(spin_system,'E')`.
- Lines 61: computes `coil` using `coil=state(spin_system,'L+','E')`.
- Lines 64-65: computes `waveform` using `waveform=linspace(parameters.sweep(1), parameters.sweep(2),parameters.nsteps)`.
- Lines 66: computes `b_axis` using `b_axis=waveform+spin_system.inter.magnet`.
- Lines 72: computes `spectrum` using `spectrum=zeros([parameters.nsteps 1],'like',1i)`.
- Lines 76: computes `spectrum(n)` using `spectrum(n)=coil'*rho`.

### Local helper functions

- Line 83: `grumble()` — `function grumble(parameters)`.
  - Representative operation: `if ~isfield(parameters,'mw_pwr')`.
  - Representative operation: `error('microwave power must be specified in parameters.mw_pwr field.')`.

## Parameters / inputs

- parameters.mw_pwr -microwave power, rad/s
- parameters.sweep -magnetic field sweep extents,
- arouind the centre field specified
- in sys.magnet, a two-element vector
- in Tesla
- parameters.nsteps -number of steps in the magnetic
- field sweep
- parameters.timestep -duration of each time step, seconds

## Outputs

- b_axis -magnetic field axis, Tesla
- spectrum -L+ observable amplitude at each
- magnetic field
- Note: this experiment should be called directly without a context.

## Implementation structure

- Time-domain rapid field scan ESR experiment, Eatons style. Syntax:
- [b_axis,spectrum]=rapidscan(spin_system,parameters)
- parameters.mw_pwr -microwave power, rad/s
- parameters.sweep -magnetic field sweep extents,
- arouind the centre field specified
- in sys.magnet, a two-element vector
- in Tesla
- parameters.nsteps -number of steps in the magnetic
- field sweep
- parameters.timestep -duration of each time step, seconds
- b_axis -magnetic field axis, Tesla
- spectrum -L+ observable amplitude at each

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `grumble()`, `operator()`, `equilibrium()`, `hamiltonian()`, `assume()`, `relaxation()`, `carrier()`, `rotframe()`, `state()`, `spectrum()`, `step()`, `waveform()`, `isfield()`, `isscalar()`.
