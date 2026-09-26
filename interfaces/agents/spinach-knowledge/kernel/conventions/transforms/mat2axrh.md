# kernel/conventions/transforms/mat2axrh.m

- Signature: `[iso,ax,rh,eigvals]=mat2axrh(M)`

## Purpose

Computes axiality and rhombicity of a symmetric 3x3 interaction tensor from the corresponding matrix. Syntax: [iso,ax,rh]=mat2axrh(M)

## Physical / mathematical content

- Convention and tensor-transform utilities. They convert among tensor parameterisations, coordinate systems, and unit systems; the underlying mathematics is linear algebra on rank-2 tensors and rotation representations.

## Numerical / algorithmic content

- An eigenvalue problem is solved or analysed, so the file is extracting spectra, stationary states, avoided crossings, or modal structure from the effective Hamiltonian or superoperator.

## Parameters / inputs

- M -a real symmetric 3x3 matrix

## Outputs

- iso -isotropic part of the interaction, defined as
- (xx+yy+zz)/3 in terms of eigenvaues
- ax -interaction axiality, defined as 2*zz-(xx+yy)
- in terms of eigenvalues
- rh -interaction rhombicity, defined as (yy-xx) in
- terms of eigenvalues
- eigvals -interaction tensor eigenvalues in Mehring order
- Note: eigenvalues [xx yy zz] are sorted in Mehring order, that
- is xx<=yy<=zz
- Note: Euler angles are not returned because the transformation
- in question is ill-defined

## Implementation structure

- Computes axiality and rhombicity of a symmetric 3x3 interaction
- tensor from the corresponding matrix. Syntax:
- [iso,ax,rh]=mat2axrh(M)
- M -a real symmetric 3x3 matrix
- iso -isotropic part of the interaction, defined as
- (xx+yy+zz)/3 in terms of eigenvaues
- ax -interaction axiality, defined as 2*zz-(xx+yy)
- in terms of eigenvalues
- rh -interaction rhombicity, defined as (yy-xx) in
- terms of eigenvalues
- eigvals -interaction tensor eigenvalues in Mehring order
- Note: eigenvalues [xx yy zz] are sorted in Mehring order, that
