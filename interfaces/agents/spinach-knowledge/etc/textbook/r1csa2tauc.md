# etc/textbook/r1csa2tauc.m

- Signature: `tauc=r1csa2tauc(R1,del_sq,B0,isotope)`

## Purpose

Estimates the rotational correlation time from the longitudinal CSA relaxation rate. Syntax: tauc=r1csa2tauc(R1,del_sq,B0,isotope)

## Physical / mathematical content

- Chemical-shift anisotropy is present: shielding is treated as a second-rank tensor whose orientation relative to the field or rotor axis modulates line shapes and transfer dynamics.

## Numerical / algorithmic content

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
