# etc/molecules/lactate.m

- Signature: `[sys,inter]=lactate(spins)`

## Purpose

Spin system of 13C-labelled lactate with the OH protons assumed to be in rapid exchange with water. Syntax: [sys,inter,bas]=lactate(spins)

## Physical / mathematical content

## Numerical / algorithmic content

## Parameters / inputs

- spins -a cell array containing the isotopes to
- import, e.g. {'1H','13C'} or a list of
- atom labels, e.g. {'CO','CA','HA'}; the
- lists can be mixed, e.g. {'1H','CA'}

## Outputs

- sys -Spinach spin system description structure
- inter -Spinach interaction description structure

## Implementation structure

- Spin system of 13C-labelled lactate with the OH protons
- assumed to be in rapid exchange with water. Syntax:
- [sys,inter,bas]=lactate(spins)
- spins -a cell array containing the isotopes to
- import, e.g. {'1H','13C'} or a list of
- atom labels, e.g. {'CO','CA','HA'}; the
- lists can be mixed, e.g. {'1H','CA'}
- sys -Spinach spin system description structure
- inter -Spinach interaction description structure
- Check consistency
- Isotopes
- Chemical shifts (approximate)
