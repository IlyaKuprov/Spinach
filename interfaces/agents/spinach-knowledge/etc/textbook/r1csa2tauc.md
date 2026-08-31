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

## Code-derived implementation details

### Comment-guided execution stages

- Lines 27-28: Check consistency; implemented by `grumble(R1,del_sq,B0,isotope)`.
- Lines 30-31: Get the Zeeman frequency; implemented by `omega=-B0*spin(isotope)`.
- Lines 33-34: Solve the quadratic equation; implemented by `tauc(1)=(del_sq*omega^2-sqrt(del_sq^2*omega^4-225*omega^2*R1^2))/(15*omega^2*R1)`.
- Lines 37-38: Check for imaginary components; implemented by `if (~isreal(tauc(1)))&&(~isreal(tauc(2)))`.

### Control flow inferred from the code

- Line 38: conditional branch on `(~isreal(tauc(1)))&&(~isreal(tauc(2)))`.

### Key state/data transformations

- Lines 31: computes `omega` using `omega=-B0*spin(isotope)`.
- Lines 34: computes `tauc(1)` using `tauc(1)=(del_sq*omega^2-sqrt(del_sq^2*omega^4-225*omega^2*R1^2))/(15*omega^2*R1)`.
- Lines 35: computes `tauc(2)` using `tauc(2)=(del_sq*omega^2+sqrt(del_sq^2*omega^4-225*omega^2*R1^2))/(15*omega^2*R1)`.

### Local helper functions

- Line 45: `grumble()` — `function grumble(R1,del_sq,B0,isotope)`.
  - Representative operation: `if (~isnumeric(R1))||(~isreal(R1))|| (~isscalar(R1))||(R1<=0)`.
  - Representative operation: `(~isscalar(R1))||(R1<=0)`.

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
- Solve the quadratic equation

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `grumble()`, `spin()`, `tauc()`, `isscalar()`, `ischar()`.
