# kernel/utilities/idxof.m

- Signature: `idx=idxof(sys,label)`

## Purpose

Allows interaction specification by spin label rather than number. Syntax: idx=idxof(sys,label)

## Physical / mathematical content

- General mathematical and infrastructure utilities. This area contains finite differences, perturbation theory, graph algorithms, spectral densities, tensor algebra, hash/report helpers, and other reusable numerical components.

## Numerical / algorithmic content

## Parameters / inputs

- sys -Spinach input structure that
- includes a sys.labels field
- with unique labels
- label -label whose index is to be returned

## Outputs

- idx -the index of the spin, an integer

## Implementation structure

- Allows interaction specification by spin label
- rather than number. Syntax:
- idx=idxof(sys,label)
- sys -Spinach input structure that
- includes a sys.labels field
- with unique labels
- label -label whose index is to be returned
- idx -the index of the spin, an integer
- Check consistency
- Locate the label
- Check the output
- Consistency enforcement
