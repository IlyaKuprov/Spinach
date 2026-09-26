# kernel/utilities/transfermat.m

- Signature: `T=transfermat(amp_inps,amp_outs)`

## Purpose

Transfer matrix calculation for linear filters. Syntax: T=transfermat(amp_inps,amp_outs)

## Physical / mathematical content

- General mathematical and infrastructure utilities. This area contains finite differences, perturbation theory, graph algorithms, spectral densities, tensor algebra, hash/report helpers, and other reusable numerical components.

## Numerical / algorithmic content

## Parameters / inputs

- amp_inps -a matrix with amplifier input vectors as columns
- amp_outs -a matrix with amplifier output vectors as columns

## Outputs

- T -the transfer matrix, such that amp_outs=T*amp_inps
- in the least squares sense
- Note: the number of input-output vector pairs should be bigger than
- the number of elements in those vectors.

## Implementation structure

- Transfer matrix calculation for linear filters. Syntax:
- T=transfermat(amp_inps,amp_outs)
- amp_inps -a matrix with amplifier input vectors as columns
- amp_outs -a matrix with amplifier output vectors as columns
- T -the transfer matrix, such that amp_outs=T*amp_inps
- in the least squares sense
- Note: the number of input-output vector pairs should be bigger than
- the number of elements in those vectors.
- Check consistency
- Run the SVD pseudoinverse
- Consistency enforcement
- A good friend will always stab you in the front.
