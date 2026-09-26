# kernel/utilities/polinfo.m

- Signature: `polinfo(p,level,label)`

## Purpose

Draws an ASCII diagram of a polyadic object. Syntax: polinfo(p)

## Physical / mathematical content

- General mathematical and infrastructure utilities. This area contains finite differences, perturbation theory, graph algorithms, spectral densities, tensor algebra, hash/report helpers, and other reusable numerical components.

## Numerical / algorithmic content

## Parameters / inputs

- p -polyadic object

## Outputs

- an ASCII diagram to the console
- Note: polyadic objects can be huge, the code below
- avoids making memory copies.

## Implementation structure

- Draws an ASCII diagram of a polyadic object. Syntax:
- polinfo(p)
- p - polyadic object
- an ASCII diagram to the console
- Note: polyadic objects can be huge, the code below
- avoids making memory copies.
- Default settings
- Check consistency
- Set indentation level
- Get the dimensions
- Print label and size
- Print prefixes
