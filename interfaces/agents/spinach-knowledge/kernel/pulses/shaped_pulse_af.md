# kernel/pulses/shaped_pulse_af.m

- Signature: `[rho,traj,P]=shaped_pulse_af(spin_system,L0,Lx,Ly,rho,rf_frq_list,...`

## Purpose

Shaped pulse in amplitude-frequency coordinates using Fokker-Planck formalism (Eqn. 33 in http://dx.doi.org/10.1016/j.jmr.2016.07.005).

## Physical / mathematical content

- Pulse and waveform utilities. These files encode shaped RF pulses, gradient events, rotating-frame transformations, resonator response, and Lie-group integration of time-dependent driven dynamics.
- The file uses a Fokker-Planck-style enlarged state space in which spatial or orientational coordinates are promoted to extra dimensions and coupled to spin dynamics through differential operators.
- Propagation is accelerated with a Krylov-subspace method, replacing direct matrix exponentiation by projection into a much smaller Arnoldi/Lanczos-type subspace.
- Signal processing is central here: the code moves between time and frequency domains, typically using FFT conventions, apodisation, zero filling, or heterodyne frequency shifts.

## Numerical / algorithmic content

- Time propagation is explicit. In Spinach this usually means repeated application of matrix exponentials or propagator factorizations to density operators or state vectors in Hilbert/Liouville/Fokker-Planck space.
- The implementation explicitly addresses performance engineering through parallel or GPU execution, which matters because Spinach operators can become extremely large after basis expansion or powder/spatial lifting.
- A Krylov-subspace or Arnoldi construction is used to avoid forming or exponentiating very large dense propagators directly.

## Syntax

```matlab
[rho,traj,P]=shaped_pulse_af(spin_system,L0,Lx,Ly,rho,rf_frq_list,...
rf_amp_list,rf_dur_list,rf_phi,...
max_rank,method)
```

## Parameters / inputs

- L0 -drift Liouvillian that continues
- running in the background
- Lx -X projection of the RF operator
- Ly -Y projection of the RF operator
- rho -initial state vector or a horizontal
- stack thereof
- rf_frq_list -a vector of RF frequencies at each
- time slice (relative to the offsets
- and/or rotating frames that were
- used to make the background L0 that
- you have supplied, Hz
- rf_amp_list -a vector of RF amplitudes at each
- time slice, rad/s
- rf_dur_list -a vector of time slice durations,
- in seconds
- rf_phi -RF phase of the first pulse slice
- max_rank -maximum rank of the Fokker-Planck
- theory, increase until the answer
- stops changing, 2 is a good start
- method -propagation method, 'expv' for Krylov
- propagation, 'expm' for exponential
- propagation, 'evolution' for Spinach
- evolution function

## Outputs

- rho -final state vector or a stack thereof
- traj -system trajectory as a [1 x (nsteps+1)]
- cell array, the first point is the ini-
- tial condition
- P -effective pulse propagator (expensive),
- only available for the 'expm' method
- Note: the pulse is assumed to be piecewise-constant and should be
- supplied with sufficiently fine time discretisation to pro-
- perly reproduce the waveform.
- Note: make it dead certain that your freqiency has the correct
- sign; wrong sign means that the pulse hits very far away
- from your intended location. This is the principal source
- of bugs when using this function.

## Implementation structure

- Shaped pulse in amplitude-frequency coordinates using Fokker-Planck
- formalism (Eqn. 33 in http://dx.doi.org/10.1016/j.jmr.2016.07.005).
- [rho,traj,P]=shaped_pulse_af(spin_system,L0,Lx,Ly,rho,rf_frq_list,...
- rf_amp_list,rf_dur_list,rf_phi,...
- max_rank,method)
- L0 -drift Liouvillian that continues
- running in the background
- Lx -X projection of the RF operator
- Ly -Y projection of the RF operator
- rho -initial state vector or a horizontal
- stack thereof
- rf_frq_list -a vector of RF frequencies at each
