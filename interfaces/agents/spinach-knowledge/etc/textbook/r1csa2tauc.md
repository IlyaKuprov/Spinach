# etc/textbook/r1csa2tauc.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/etc/textbook/r1csa2tauc.m`
- Signature: `tauc=r1csa2tauc(R1,del_sq,B0,isotope)`
- Total lines: 67

## Purpose

Estimates the rotational correlation time from the longitudinal CSA relaxation rate. Syntax: tauc=r1csa2tauc(R1,del_sq,B0,isotope)

## Physical / mathematical content

- This file belongs to the `etc` part of Spinach. Its role should be read together with nearby files in the same directory, which usually share a common physical regime or infrastructure purpose.
- Chemical-shift anisotropy is present: shielding is treated as a second-rank tensor whose orientation relative to the field or rotor axis modulates line shapes and transfer dynamics.

## Numerical / algorithmic content

- The file contains an explicit `grumble(...)` validator, which is Spinach convention for front-loading dimension, type, and regime checks before expensive linear-algebra work begins.
- The file also defines local helper function(s): `grumble()`. This usually means the public entry point is supported by tightly coupled validation or helper logic kept private to the file.

## Parameters / inputs

- R1 -longitudinal relaxation rate, Hz
- del_sq -second rank invariant of the CSA,
- see blinv.m function
- B0 -magnetic field, Tesla
- isotope -isotope specification string, e.g. '1H'

## Outputs

- tauc -rotational correlation time, seconds

## Implementation structure

- Estimates the rotational correlation time from the
- longitudinal CSA relaxation rate. Syntax:
- tauc=r1csa2tauc(R1,del_sq,B0,isotope)
- R1 -longitudinal relaxation rate, Hz
- del_sq -second rank invariant of the CSA,
- see blinv.m function
- B0 -magnetic field, Tesla
- isotope -isotope specification string, e.g. '1H'
- tauc -rotational correlation time, seconds
- Check consistency
- Get the Zeeman frequency
- Solve the quadratic equation, larger root first to avoid cancellation

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `grumble()`, `spin()`, `tauc()`, `isscalar()`, `ischar()`.
