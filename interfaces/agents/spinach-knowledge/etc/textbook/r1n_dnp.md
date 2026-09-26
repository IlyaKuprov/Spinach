# etc/textbook/r1n_dnp.m

- Signature: `R1n=r1n_dnp(B0,T,g,T1e,T1n_bulk,r,bet)`

## Purpose

A simple model for nuclear longitudnal relaxation rate dies to the presence of an unpaired electron in cryogenic DNP settings. [literature reference goes here]

## Physical / mathematical content

## Numerical / algorithmic content

## Syntax

```matlab
R1n=r1n_dnp(B0,T,g,T1e,T1n_bulk,r,bet)
```

## Parameters / inputs

- B0 -main magnet field, Tesla
- T -absolute temperature, Kelvin
- g -electron g-factor, Bohr
- magneton units
- T1e -electron longitudinal rela-
- xation time, seconds
- T1n_bulk -nuclear longitudinal relaxa-
- tion time far away from the
- electron
- r -electron-nuclear distance,
- Angstrom
- bet -angle between the magnet field
- and the electron-nuclear direc-
- tion, radians

## Outputs

- R1n -nuclear relaxation rate, Hz

## Implementation structure

- A simple model for nuclear longitudnal relaxation
- rate dies to the presence of an unpaired electron
- in cryogenic DNP settings.
- [literature reference goes here]
- R1n=r1n_dnp(B0,T,g,T1e,T1n_bulk,r,bet)
- B0 -main magnet field, Tesla
- T -absolute temperature, Kelvin
- g -electron g-factor, Bohr
- magneton units
- T1e -electron longitudinal rela-
- xation time, seconds
- T1n_bulk -nuclear longitudinal relaxa-
