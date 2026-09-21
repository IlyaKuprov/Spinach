# kernel/pulses/grad_pulse.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/kernel/pulses/grad_pulse.m`
- Signature: `rho=grad_pulse(spin_system,L,rho,g_amp,s_len,g_dur,s_fac)`
- Total lines: 114

## Purpose

Emulates the effect of a gradient pulse on the sample average density matrix using Edwards formalism. It is assumed that the effect of dif- fusion is negligible, that the gradient is linear, and that it is an- tisymmetric about the middle of the sample. Syntax: rho=grad_pulse(spin_system,rho,g_amp,s_len,g_dur,s_fac)

## Physical / mathematical content

- Pulse and waveform utilities. These files encode shaped RF pulses, gradient events, rotating-frame transformations, resonator response, and Lie-group integration of time-dependent driven dynamics.

## Numerical / algorithmic content

- Time propagation is explicit. In Spinach this usually means repeated application of matrix exponentials or propagator factorizations to density operators or state vectors in Hilbert/Liouville/Fokker-Planck space.
- The file contains an explicit `grumble(...)` validator, which is Spinach convention for front-loading dimension, type, and regime checks before expensive linear-algebra work begins.
- The file also defines local helper function(s): `grumble()`. This usually means the public entry point is supported by tightly coupled validation or helper logic kept private to the file.

## Parameters / inputs

- rho -spin system state vector
- L -system Liouvillian
- g_amp -gradient amplitude, Gauss/cm
- s_len -sample length, cm
- g_dur -gradient pulse duration, seconds
- s_fac -gradient shape factor, use 1 for
- square gradient pulses

## Outputs

- rho -spin system state vector, integrated over
- the spatial coordinate
- Note: the function integrates over sample coordinates -subsequent gra-
- dient pulses would not refocus the magnetization that it has de-
- focused. To simulate a gradient sandwich, use grad_sandw.m func-
- tion. More information on the subject is available in Luke's pa-
- per (http://dx.doi.org/10.1016/j.jmr.2014.01.011).
- Note: this function is OK for standalone crusher gradients; for more
- sophisticated gradient work, use the imaging context.

## Implementation structure

- Emulates the effect of a gradient pulse on the sample average density
- matrix using Edwards formalism. It is assumed that the effect of dif-
- fusion is negligible, that the gradient is linear, and that it is an-
- tisymmetric about the middle of the sample. Syntax:
- rho=grad_pulse(spin_system,rho,g_amp,s_len,g_dur,s_fac)
- rho -spin system state vector
- L -system Liouvillian
- g_amp -gradient amplitude, Gauss/cm
- s_len -sample length, cm
- g_dur -gradient pulse duration, seconds
- s_fac -gradient shape factor, use 1 for
- square gradient pulses

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `grumble()`, `report()`, `carrier()`, `cheap_norm()`, `evolution()`, `speye()`, `aux_rho()`, `ismember()`, `isscalar()`.
