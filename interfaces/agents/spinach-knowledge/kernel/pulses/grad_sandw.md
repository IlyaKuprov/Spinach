# kernel/pulses/grad_sandw.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/kernel/pulses/grad_sandw.m`
- Signature: `rho=grad_sandw(spin_system,L,rho,P,g_amps,s_len,g_durs,s_facs)`
- Total lines: 113

## Purpose

Emulates the effect of a gradient sandwich on the sample average density matrix using Edwards formalism. It is assumed that the effect of diffusi- on is negligible, that the gradients are linear, and that they are anti- symmetric about the middle of the sample. Syntax: rho=grad_sandw(spin_system,L,rho,P,g_amps,s_len,g_durs,s_facs)

## Physical / mathematical content

- Pulse and waveform utilities. These files encode shaped RF pulses, gradient events, rotating-frame transformations, resonator response, and Lie-group integration of time-dependent driven dynamics.

## Numerical / algorithmic content

- Time propagation is explicit. In Spinach this usually means repeated application of matrix exponentials or propagator factorizations to density operators or state vectors in Hilbert/Liouville/Fokker-Planck space.
- The file contains an explicit `grumble(...)` validator, which is Spinach convention for front-loading dimension, type, and regime checks before expensive linear-algebra work begins.
- The file also defines local helper function(s): `grumble()`. This usually means the public entry point is supported by tightly coupled validation or helper logic kept private to the file.

## Code-derived implementation details

### Comment-guided execution stages

- Lines 48-49: Check consistency; implemented by `grumble(spin_system,L,rho,P,g_amps,s_len,g_durs,s_facs)`.
- Lines 51-52: Inform the user; implemented by `report(spin_system,'computing the effect of a gradient sandwich ')`.
- Lines 54-55: Get effective gradient operators (WARNING: shifts are ignored); implemented by `R=carrier(spin_system,'all')/spin_system.inter.magnet`.
- Lines 59-61: Check approximation applicability; implemented by `if (cheap_norm(L*G1-G1*L)>1e-6)|| (cheap_norm(L*G2-G2*L)>1e-6)`.
- Lines 65-66: Run background Liouvillian evolution during first gradient; implemented by `rho=evolution(spin_system,L,[],rho,g_durs(1),1,'final')`.
- Lines 68-69: Account for rescaling of sample integral into [0,1]; implemented by `rho=evolution(spin_system,G1,[],rho,-1/2,1,'final')`.
- Lines 71-72: Evolution under the gradient pulses and P; implemented by `aux_mat=[-G2, 1i*P; 0*G1, G1]; aux_rho=[0*rho; rho]`.
- Lines 75-76: Map back into the state vector space; implemented by `rho=aux_rho(1:end/2,:)`.
- Lines 78-79: Undo the evolution overshoot; implemented by `rho=evolution(spin_system,G2,[],rho,1/2,1,'final')`.
- Lines 81-82: Run background Liouvillian evolution during second gradient; implemented by `rho=evolution(spin_system,L,[],rho,g_durs(2),1,'final')`.

### Control flow inferred from the code

- Line 60: conditional branch on `(cheap_norm(L*G1-G1*L)>1e-6)||`.

### Key state/data transformations

- Lines 55: computes `R` using `R=carrier(spin_system,'all')/spin_system.inter.magnet`.
- Lines 56: computes `G1` using `G1=1e-4*s_facs(1)*g_amps(1)*s_len*g_durs(1)*R`.
- Lines 57: computes `G2` using `G2=1e-4*s_facs(2)*g_amps(2)*s_len*g_durs(2)*R`.
- Lines 66: computes `rho` using `rho=evolution(spin_system,L,[],rho,g_durs(1),1,'final')`.
- Lines 72: computes `aux_mat` using `aux_mat=[-G2, 1i*P; 0*G1, G1]; aux_rho=[0*rho; rho]`.
- Lines 73: computes `aux_rho` using `aux_rho=evolution(spin_system,aux_mat,[],aux_rho,1,1,'final')`.

### Local helper functions

- Line 87: `grumble()` — `function grumble(spin_system,L,rho,P,g_amps,s_len,g_durs,s_facs)`.
  - Representative operation: `if ~ismember(spin_system.bas.formalism,{'sphten-liouv','zeeman-liouv'})`.
  - Representative operation: `error('this function is only applicable in Liouville space.')`.

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

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `grumble()`, `report()`, `carrier()`, `s_facs()`, `g_amps()`, `g_durs()`, `cheap_norm()`, `evolution()`, `aux_rho()`, `ismember()`, `isscalar()`, `any()`.
