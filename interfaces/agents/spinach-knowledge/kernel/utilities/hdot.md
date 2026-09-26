# kernel/utilities/hdot.m

- Signature: `H=hdot(A,B)`

## Purpose

Hadamard route to Frobenius matrix product. Useful as a replacement for trace(A'*B) because trace(A'*B)=hadm(conj(A),B) and the latter only needs O(n^2) multiplications as com- pared to O(n^3) for trace(A'*B). Syntax: H=hdot(A,B)

## Physical / mathematical content

- General mathematical and infrastructure utilities. This area contains finite differences, perturbation theory, graph algorithms, spectral densities, tensor algebra, hash/report helpers, and other reusable numerical components.

## Numerical / algorithmic content

## Parameters / inputs

- A,B -square matrices of the same size

## Outputs

- H -Frobenius inner product of A and B

## Implementation structure

- Hadamard route to Frobenius matrix product. Useful as a
- replacement for trace(A'*B) because
- trace(A'*B)=hadm(conj(A),B)
- and the latter only needs O(n^2) multiplications as com-
- pared to O(n^3) for trace(A'*B). Syntax:
- H=hdot(A,B)
- A,B -square matrices of the same size
- H -Frobenius inner product of A and B
- Check consistency
- Do the calculation
- Consistency enforcement
- An infinite number of mathematicians walk into a bar. The first one
