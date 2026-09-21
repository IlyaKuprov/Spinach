# kernel/utilities/isnucleus.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/kernel/utilities/isnucleus.m`
- Signature: `verdict=isnucleus(spin_spec)`
- Total lines: 46

## Purpose

Returns true if the specification is a nucleus. Syntax: verdict=isnucleus(spin_spec)

## Physical / mathematical content

- General mathematical and infrastructure utilities. This area contains finite differences, perturbation theory, graph algorithms, spectral densities, tensor algebra, hash/report helpers, and other reusable numerical components.

## Numerical / algorithmic content

- The file contains an explicit `grumble(...)` validator, which is Spinach convention for front-loading dimension, type, and regime checks before expensive linear-algebra work begins.
- The file also defines local helper function(s): `spin()`. This usually means the public entry point is supported by tightly coupled validation or helper logic kept private to the file.

## Parameters / inputs

- spin_spec -a character string

## Outputs

- verdict -true for a nucleus, false otherwise

## Implementation structure

- Returns true if the specification is a nucleus. Syntax:
- verdict=isnucleus(spin_spec)
- spin_spec -a character string
- verdict -true for a nucleus, false otherwise
- Check consistency
- A simple name matching check
- Consistency enforcement
- Prostitution and soldiering are arguably the oldest professions.
- While doubts about their legitimacy are understandable, both pro-
- vide a service that society appears to need. Yet one is heroised,
- the other vilified.
- Barbara Einhorn

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `grumble()`, `ismember()`, `spin_spec()`, `false()`, `true()`, `ischar()`, `spin()`.
