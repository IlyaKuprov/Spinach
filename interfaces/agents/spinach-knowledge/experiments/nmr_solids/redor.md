# experiments/nmr_solids/redor.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/experiments/nmr_solids/redor.m`
- Signature: `redor_curve=redor(spin_system,parameters,H,R,K)`
- Total lines: 251

## Purpose

Rotational-echo double-resonance (REDOR) experiment with ideal hard pi pulses. The observed channel is refocused once per rotor period, and the dephasing channel is pulsed at rotor-period boundaries. The sequence reports the full rotational echo, the dephased echo, and their difference as a function of the number of rotor cycles. To be called from singlerot context. Further information in:

## Physical / mathematical content

- Solid-state pulse sequence implementations. The core ingredients are anisotropic Hamiltonians, rotor synchronisation, cross-polarisation, recoupling/decoupling, and powder or rotor-stack propagation.

## Numerical / algorithmic content

- Time propagation is explicit. In Spinach this usually means repeated application of matrix exponentials or propagator factorizations to density operators or state vectors in Hilbert/Liouville/Fokker-Planck space.
- The file contains an explicit `grumble(...)` validator, which is Spinach convention for front-loading dimension, type, and regime checks before expensive linear-algebra work begins.
- The file also defines local helper function(s): `grumble()`. This usually means the public entry point is supported by tightly coupled validation or helper logic kept private to the file.

## Code-derived implementation details

### Comment-guided execution stages

- Lines 70-71: Set default decoupling; implemented by `if ~isfield(parameters,'decouple')`.
- Lines 75-76: Set default refocusing phase; implemented by `if ~isfield(parameters,'refocus_phase')`.
- Lines 80-81: Set default dephasing pulse phases; implemented by `if ~isfield(parameters,'pulse_phase')`.
- Lines 85-86: Check consistency; implemented by `grumble(spin_system,parameters,H,R,K)`.
- Lines 88-89: Compose the Liouvillian; implemented by `L=H+1i*R+1i*K`.
- Lines 91-92: Apply analytical decoupling if requested; implemented by `[L,rho0]=decouple(spin_system,L,parameters.rho0,parameters.decouple)`.
- Lines 94-95: Get observed-channel refocusing operator; implemented by `Ip=operator(spin_system,'L+',parameters.spins{1})`.
- Lines 101-102: Get dephasing-channel pulse operators; implemented by `Ip=operator(spin_system,'L+',parameters.spins{2})`.
- Lines 107-108: Precompute phased dephasing-channel pulse operators; implemented by `deph_opers=cell(size(parameters.pulse_phase))`.
- Lines 114-115: Get rotor-synchronised timing; implemented by `rotor_period=1/abs(parameters.rate)`.
- Lines 118-119: Preallocate echo arrays; implemented by `s0_echo=zeros(1,max(parameters.ncycles)+1)`.
- Lines 122-123: Store the zero-time point; implemented by `s0_echo(1)=parameters.coil'*rho0`.
- Lines 126-127: Initialise reference and dephased trajectories; implemented by `rho_s0=rho0; rho_s=rho0; pulse_count=0`.
- Lines 129-130: Step through the requested REDOR evolution window; implemented by `for n=1:max(parameters.ncycles)`.
- Lines 132-133: Propagate both echoes over the first half-period; implemented by `rho_s0=step(spin_system,L,rho_s0,half_period)`.
- Lines 136-137: Refocus the observed channel; implemented by `rho_s0=step(spin_system,ref_oper,rho_s0,pi)`.
- Lines 140-141: Propagate both echoes over the second half-period; implemented by `rho_s0=step(spin_system,L,rho_s0,half_period)`.
- Lines 144-145: Apply the dephasing-channel pi pulse to the S echo; implemented by `pulse_count=pulse_count+1`.

### Control flow inferred from the code

- Line 71: conditional branch on `~isfield(parameters,'decouple')`.
- Line 76: conditional branch on `~isfield(parameters,'refocus_phase')`.
- Line 81: conditional branch on `~isfield(parameters,'pulse_phase')`.
- Line 109: `for` loop over `n=1:numel(parameters.pulse_phase)`.
- Line 130: `for` loop over `n=1:max(parameters.ncycles)`.

### Key state/data transformations

- Lines 72: computes `parameters.decouple` using `parameters.decouple={}`.
- Lines 77: computes `parameters.refocus_phase` using `parameters.refocus_phase=0`.
- Lines 82: computes `parameters.pulse_phase` using `parameters.pulse_phase=[0 pi/2]`.
- Lines 89: computes `L` using `L=H+1i*R+1i*K`.
- Lines 92: computes `[L,rho0]` using `[L,rho0]=decouple(spin_system,L,parameters.rho0,parameters.decouple)`.
- Lines 95: computes `Ip` using `Ip=operator(spin_system,'L+',parameters.spins{1})`.
- Lines 96: computes `Ix` using `Ix=(Ip+Ip')/2; Iy=(Ip-Ip')/2i`.
- Lines 97-98: computes `ref_oper` using `ref_oper=cos(parameters.refocus_phase)*Ix+ sin(parameters.refocus_phase)*Iy`.
- Lines 105: computes `Iy` using `Iy=kron(speye(parameters.spc_dim),Iy)`.
- Lines 108: computes `deph_opers` using `deph_opers=cell(size(parameters.pulse_phase))`.
- Lines 110-111: computes `deph_opers{n}` using `deph_opers{n}=cos(parameters.pulse_phase(n))*Ix+ sin(parameters.pulse_phase(n))*Iy`.
- Lines 115: computes `rotor_period` using `rotor_period=1/abs(parameters.rate)`.
- Lines 116: computes `half_period` using `half_period=rotor_period/2`.
- Lines 119: computes `s0_echo` using `s0_echo=zeros(1,max(parameters.ncycles)+1)`.
- Lines 120: computes `s_echo` using `s_echo=zeros(1,max(parameters.ncycles)+1)`.
- Lines 123: computes `s0_echo(1)` using `s0_echo(1)=parameters.coil'*rho0`.
- Lines 124: computes `s_echo(1)` using `s_echo(1)=parameters.coil'*rho0`.
- Lines 127: computes `rho_s0` using `rho_s0=rho0; rho_s=rho0; pulse_count=0`.

### Local helper functions

- Line 165: `grumble()` — `function grumble(spin_system,parameters,H,R,K)`.
  - Representative operation: `if ~ismember(spin_system.bas.formalism,{'sphten-liouv'})`.
  - Representative operation: `error('this function is only available for sphten-liouv formalism.')`.

## Syntax

```matlab
redor_curve=redor(spin_system,parameters,H,R,K)
```

## Parameters / inputs

- parameters.spins -observed and dephasing spins,
- e.g. {'13C','15N'} for 13C{15N}
- REDOR
- parameters.ncycles -row vector with numbers of rotor
- cycles in the REDOR evolution time
- parameters.rate -MAS rate in Hz
- parameters.rho0 -initial state vector, usually
- transverse magnetisation on the
- observed spin
- parameters.coil -detection state vector, usually on
- the observed spin
- parameters.decouple -nuclei to decouple analytically
- during REDOR evolution, supplied as
- a cell array of isotope strings
- parameters.refocus_phase -phase of the observed-channel pi
- refocusing pulses, radians; defaults
- to 0
- parameters.pulse_phase -phases of the dephasing-channel pi
- pulses, radians; the vector is cycled
- through the pulse train, and defaults
- to [0 pi/2]
- H -Hamiltonian matrix, received from context function
- R -relaxation superoperator, received from context function
- K -kinetics superoperator, received from context function

## Outputs

- redor_curve(1,:) -full echo S0, with observed-channel
- refocusing pulses only
- redor_curve(2,:) -dephased echo S, with observed-channel
- refocusing and dephasing-channel pi
- pulses
- redor_curve(3,:) -REDOR difference, S0-S

## Implementation structure

- Rotational-echo double-resonance (REDOR) experiment with ideal
- hard pi pulses. The observed channel is refocused once per rotor
- period, and the dephasing channel is pulsed at rotor-period
- boundaries. The sequence reports the full rotational echo, the
- dephased echo, and their difference as a function of the number of
- rotor cycles. To be called from singlerot context. Further
- information in:
- redor_curve=redor(spin_system,parameters,H,R,K)
- parameters.spins -observed and dephasing spins,
- e.g. {'13C','15N'} for 13C{15N}
- REDOR
- parameters.ncycles -row vector with numbers of rotor

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `isfield()`, `grumble()`, `decouple()`, `operator()`, `speye()`, `s0_echo()`, `s_echo()`, `step()`, `ismember()`, `ismatrix()`, `all()`, `iscell()`, `any()`, `cellfun()`, `isrow()`, `isscalar()`.
