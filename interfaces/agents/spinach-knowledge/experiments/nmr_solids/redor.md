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
