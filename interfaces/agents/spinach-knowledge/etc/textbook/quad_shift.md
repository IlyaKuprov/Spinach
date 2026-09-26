# etc/textbook/quad_shift.m

- Signature: `delta=quad_shift(Cq,eta,v0,S,m)`

## Purpose

Second order shift of the centre of gravity of the powder pattern of |S,m> to |S,m-1> transition in the NMR spectrum of a quadrupo- lar nucleus with spin S. Equation (3) from

## Physical / mathematical content

- Quadrupolar physics is relevant: nuclei with spin > 1/2 interact with the electric field gradient tensor, introducing second-rank anisotropy, asymmetry, and overtone or MQ phenomena.

## Numerical / algorithmic content

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
