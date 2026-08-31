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

## Code-derived implementation details

### Comment-guided execution stages

- Lines 44-45: Check consistency; implemented by `grumble(spin_system,L,rho,g_amp,s_len,g_dur,s_fac)`.
- Lines 47-48: Inform the user; implemented by `report(spin_system,'computing the effect of a gradient pulse ')`.
- Lines 50-51: Get effective gradient operator (WARNING: shifts are ignored); implemented by `G=1e-4*s_fac*g_amp*s_len*g_dur*carrier(spin_system,'all')/spin_system.inter.magnet`.
- Lines 53-54: Check approximation applicability; implemented by `if cheap_norm(L*G-G*L)>1e-6`.
- Lines 58-59: Run background Liouvillian evolution; implemented by `rho=evolution(spin_system,L,[],rho,g_dur,1,'final')`.
- Lines 61-62: Account for rescaling of sample integral into [0,1]; implemented by `rho=evolution(spin_system,G,[],rho,-1/2,1,'final')`.
- Lines 64-65: Compute sample integral using auxiliary matrix method; implemented by `aux_mat=[0*G, 1i*speye(size(G)); 0*G, G]; aux_rho=[0*rho; rho]`.
- Lines 68-69: Map back into the state vector space; implemented by `rho=aux_rho(1:end/2,:)`.
- Lines 71-72: Inform the user; implemented by `report(spin_system,'gradient propagation done.')`.

### Control flow inferred from the code

- Line 54: conditional branch on `cheap_norm(L*G-G*L)>1e-6`.

### Key state/data transformations

- Lines 51: computes `G` using `G=1e-4*s_fac*g_amp*s_len*g_dur*carrier(spin_system,'all')/spin_system.inter.magnet`.
- Lines 59: computes `rho` using `rho=evolution(spin_system,L,[],rho,g_dur,1,'final')`.
- Lines 65: computes `aux_mat` using `aux_mat=[0*G, 1i*speye(size(G)); 0*G, G]; aux_rho=[0*rho; rho]`.
- Lines 66: computes `aux_rho` using `aux_rho=evolution(spin_system,aux_mat,[],aux_rho,1,1,'final')`.

### Local helper functions

- Line 77: `grumble()` — `function grumble(spin_system,L,rho,g_amp,s_len,g_dur,s_fac)`.
  - Representative operation: `if ~ismember(spin_system.bas.formalism,{'sphten-liouv','zeeman-liouv'})`.
  - Representative operation: `error('this function is only applicable in Liouville space.')`.

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
