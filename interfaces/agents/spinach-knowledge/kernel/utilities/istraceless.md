# kernel/utilities/istraceless.m

- Signature: `A=istraceless(M)`

## Purpose

A floating-point precision consistent check for whether a particular matrix is traceless. Syntax: A=istraceless(M)

## Physical / mathematical content

- General mathematical and infrastructure utilities. This area contains finite differences, perturbation theory, graph algorithms, spectral densities, tensor algebra, hash/report helpers, and other reusable numerical components.

## Numerical / algorithmic content

## Parameters / inputs

- M -a matrix of any dimension

## Outputs

- A -true if the matrix is traceless to ap-
- propriate precision, false otherwise

## Implementation structure

- A floating-point precision consistent check for whether
- a particular matrix is traceless. Syntax:
- A=istraceless(M)
- M -a matrix of any dimension
- A -true if the matrix is traceless to ap-
- propriate precision, false otherwise
- Check consistency
- Working precision
- Cheapest norm of M
- Decide if M is traceless
- Consistency enforcement
- The College asked me to chair the Size and Shape
