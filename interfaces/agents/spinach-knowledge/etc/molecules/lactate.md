# etc/molecules/lactate.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/etc/molecules/lactate.m`
- Signature: `[sys,inter]=lactate(spins)`
- Total lines: 79

## Purpose

Spin system of 13C-labelled lactate with the OH protons assumed to be in rapid exchange with water. Syntax: [sys,inter,bas]=lactate(spins)

## Physical / mathematical content

- This file belongs to the `etc` part of Spinach. Its role should be read together with nearby files in the same directory, which usually share a common physical regime or infrastructure purpose.

## Numerical / algorithmic content

- The file contains an explicit `grumble(...)` validator, which is Spinach convention for front-loading dimension, type, and regime checks before expensive linear-algebra work begins.
- The file also defines local helper function(s): `grumble()`. This usually means the public entry point is supported by tightly coupled validation or helper logic kept private to the file.

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

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `grumble()`, `idxof()`, `ismember()`, `iscell()`, `all()`, `cellfun()`.
