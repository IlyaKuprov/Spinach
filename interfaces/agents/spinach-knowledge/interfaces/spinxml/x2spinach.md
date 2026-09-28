# interfaces/spinxml/x2spinach.m

- Signature: `[sys,inter]=x2spinach(filename,shielding_refs)`

## Purpose

Reads SpinXML and builds Spinach spin-system and interaction input structures.

## Physical / mathematical content

## Numerical / algorithmic content

## Parameters / inputs

- filename — a non-empty character string naming an existing SpinXML file
- shielding_refs — isotope-specific reference values used when the file supplies absolute shieldings, e.g. `{{'1H',31.5},{'13C',189.7}}`. This is needed because Spinach uses chemical shifts; pass an empty cell array when the file supplies shifts. Isotope names must be unique, and each reference value numeric.

## Outputs

- sys — Spinach spin-system input structure, including isotope and optional label data
- inter — Spinach interaction input structure, including coordinates and parsed interaction terms
- WARNING: this function assumes the SpinXML file has passed validation against the schema at http://spindynamics.org/SpinXML.php

## Implementation structure

- Validates the filename and shielding_refs, then parses the XML with parsexml().
- Makes a first pass over spin entries to collect IDs, isotopes, optional labels, and coordinates; orders spins by ID and requires IDs to start at 1 without gaps.
- Makes a second pass to parse interaction kinds, spin pairs, scalar or tensor specifications, units, and orientation data into the Spinach interaction fields.
- Uses shielding_refs to convert absolute shielding data to the chemical-shift input Spinach requires.
