# kernel/pulses/grad_sandw.m

- Signature: `rho=grad_sandw(spin_system,L,rho,P,g_amps,s_len,g_durs,s_facs)`

## Purpose

Emulates the effect of a gradient sandwich on the sample average density matrix using Edwards formalism. It is assumed that the effect of diffusi- on is negligible, that the gradients are linear, and that they are anti- symmetric about the middle of the sample. Syntax: rho=grad_sandw(spin_system,L,rho,P,g_amps,s_len,g_durs,s_facs)

## Physical / mathematical content

- Pulse and waveform utilities. These files encode shaped RF pulses, gradient events, rotating-frame transformations, resonator response, and Lie-group integration of time-dependent driven dynamics.

## Numerical / algorithmic content

- Time propagation is explicit. In Spinach this usually means repeated application of matrix exponentials or propagator factorizations to density operators or state vectors in Hilbert/Liouville/Fokker-Planck space.

## Parameters / inputs

- rho -spin system state vector
- L -system Liouvillian
- P -total propagator for all events happening
- between the two gradients
- g_amps -row vector containing the amplitudes of
- the two gradients, Gauss/cm
- s_len -sample length, cm
- g_durs -row vector containing the durations of
- the two gradients, seconds
- s_facs -shape factors of the two gradients, use
- [1 1] for square gradient pulses

## Outputs

- rho -spin system state vector, integrated over
- the spatial coordinate
- Note: the function integrates over sample coordinates -subsequent gra-
- dient pulses would not refocus the magnetization that it has left
- defocused. More information on the subject is available in Luke's
- paper (http://dx.doi.org/10.1016/j.jmr.2014.01.011).
- Note: this function is OK for standalone gradient pairs; for more
- sophisticated gradient work, use the imaging context.

## Implementation structure

- Emulates the effect of a gradient sandwich on the sample average density
- matrix using Edwards formalism. It is assumed that the effect of diffusi-
- on is negligible, that the gradients are linear, and that they are anti-
- symmetric about the middle of the sample. Syntax:
- rho=grad_sandw(spin_system,L,rho,P,g_amps,s_len,g_durs,s_facs)
- rho -spin system state vector
- L -system Liouvillian
- P -total propagator for all events happening
- between the two gradients
- g_amps -row vector containing the amplitudes of
- the two gradients, Gauss/cm
- s_len -sample length, cm
