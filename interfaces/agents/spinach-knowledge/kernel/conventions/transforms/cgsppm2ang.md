# kernel/conventions/transforms/cgsppm2ang.m

- Signature: `ang=cgsppm2ang(cgsppm)`

## Purpose

Converts magnetic susceptibility from the cgs-ppm (aka cm^3/mol) units quoted by quantum chemistry packages into Angstrom^3 units required by Spinach pseudocontact shift functionality. Syntax: ang=cgsppm2ang(cgsppm)

## Physical / mathematical content

- Convention and tensor-transform utilities. They convert among tensor parameterisations, coordinate systems, and unit systems; the underlying mathematics is linear algebra on rank-2 tensors and rotation representations.

## Numerical / algorithmic content

## Parameters / inputs

- cgsppm -any numerical array of susceptibility
- values in cgs-ppm

## Outputs

- ang -array of the same size with suscepti-
- bility values in cubic Angstrom

## Implementation structure

- Converts magnetic susceptibility from the cgs-ppm (aka cm^3/mol) units
- quoted by quantum chemistry packages into Angstrom^3 units required by
- Spinach pseudocontact shift functionality. Syntax:
- ang=cgsppm2ang(cgsppm)
- cgsppm -any numerical array of susceptibility
- values in cgs-ppm
- ang -array of the same size with suscepti-
- bility values in cubic Angstrom
- Check consistency
- Do the calculation
- Consistency enforcement
- "What is this thing, anyway?" said the Dean, inspecting the implement in
