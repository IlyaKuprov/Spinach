# etc/textbook/quad_shift.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/etc/textbook/quad_shift.m`
- Signature: `delta=quad_shift(Cq,eta,v0,S,m)`
- Total lines: 73

## Purpose

Second order shift of the centre of gravity of the powder pattern of |S,m> to |S,m-1> transition in the NMR spectrum of a quadrupo- lar nucleus with spin S. Equation (3) from

## Physical / mathematical content

- This file belongs to the `etc` part of Spinach. Its role should be read together with nearby files in the same directory, which usually share a common physical regime or infrastructure purpose.
- Quadrupolar physics is relevant: nuclei with spin > 1/2 interact with the electric field gradient tensor, introducing second-rank anisotropy, asymmetry, and overtone or MQ phenomena.

## Numerical / algorithmic content

- The file contains an explicit `grumble(...)` validator, which is Spinach convention for front-loading dimension, type, and regime checks before expensive linear-algebra work begins.
- The file also defines local helper function(s): `grumble()`. This usually means the public entry point is supported by tightly coupled validation or helper logic kept private to the file.

## Syntax

```matlab
delta=quad_shift(Cq,eta,v0,S,m)
```

## Parameters / inputs

- Cq -quadrupolar constant, Hz
- eta -quadrupolar asymmetry parameter
- v0 -Larmor frequency of the nucleus, Hz
- S -spin quantum number of the nucleus
- m -projection quantum number of the
- starting energy level

## Outputs

- delta -quadrupolar shift in ppm
- Note: a few papers contain an incorrect version of this expressi-
- on; the one used here was tested against pure numerics and
- found to be correct.

## Implementation structure

- Second order shift of the centre of gravity of the powder pattern
- of |S,m> to |S,m-1> transition in the NMR spectrum of a quadrupo-
- lar nucleus with spin S. Equation (3) from
- delta=quad_shift(Cq,eta,v0,S,m)
- Cq -quadrupolar constant, Hz
- eta -quadrupolar asymmetry parameter
- v0 -Larmor frequency of the nucleus, Hz
- S -spin quantum number of the nucleus
- m -projection quantum number of the
- starting energy level
- delta -quadrupolar shift in ppm
- Note: a few papers contain an incorrect version of this expressi-

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `grumble()`, `isscalar()`.
