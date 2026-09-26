# experiments/pseudocon/lpcs.m

- Signature: `theo_pcs=lpcs(nxyz,mxyz,ranks,Ilm,chi)`

## Purpose

Computes PCS from the multipole moments of the paramagnetic centre probability density as described in Equation 33 is used under the assumption that all nuclei are outside the bounding sphere shown in Figure 1 of the paper.

## Physical / mathematical content

- Paramagnetic-pseudocontact inference routines. The mathematics includes inverse problems, tensor parameterisation, interpolation, and regularisation.

## Numerical / algorithmic content

## Syntax

```matlab
theo_pcs=lpcs(nxyz,mxyz,ranks,Ilm,chi)
```

## Parameters / inputs

- nxyz -nuclear coordinates as [x y z] with multiple
- rows, at which PCS is to be evaluated, in
- Angstroms.
- mxyz -paramagnetic centre coordinates as [x y z],
- in Angstroms.
- ranks -array of multipole ranks supplied in Ilm
- Ilm -{[],[]} array of numbers corresponding to
- the integrals
- Int[rho(r,theta,phi)*Y_lm(theta,phi)*r^l*d^3r]
- for L=0 Ilm=N/2/sqrt(pi)
- for L=1 Ilm=[imag(I11) I10 real(I11)]
- for L=2 Ilm=[imag(I22) imag(I21) I20 ...
- real(I21) real(I22)]
- et cetera.
- chi -the 3x3 matrix or the five independent elements
- of the susceptibility tensor in cubic Angstroms
- Output:
- theo_pcs -predicted pseudocontact shift (in ppm) at
- each of the nuclei.

## Implementation structure

- Computes PCS from the multipole moments of the paramagnetic
- centre probability density as described in
- Equation 33 is used under the assumption that all nuclei are
- outside the bounding sphere shown in Figure 1 of the paper.
- theo_pcs=lpcs(nxyz,mxyz,ranks,Ilm,chi)
- nxyz -nuclear coordinates as [x y z] with multiple
- rows, at which PCS is to be evaluated, in
- Angstroms.
- mxyz -paramagnetic centre coordinates as [x y z],
- in Angstroms.
- ranks -array of multipole ranks supplied in Ilm
- Ilm -{[],[]} array of numbers corresponding to
