# kernel/states/stateinfo.m

- Signature: `stateinfo(spin_system,rho,npops)`

## Purpose

Prints the state-vector 2-norm and the basis states associated with the `npops` largest-magnitude coefficients in `rho`, in descending order. Syntax: `stateinfo(spin_system,rho,npops)`

## Physical / mathematical content

- Requires a spherical-tensor basis; each reported label identifies a direct-product spherical-tensor component.

## Parameters / inputs

- rho -state vector
- npops -number of largest-magnitude coefficients to print

## Outputs

- This function prints a summary of the state composition to the con-
- sole in the following format:
- (L1,M1) (L2,M2) ... (Ln,Mn) coefficient number
- This corresponds to the direct product of single-spin irreducible
- spherical tensors with the specified indices, its coefficient in
- the linear combination, and the number of the corresponding state
- in the basis set.
- Note: this function requires a spherical tensor basis set.
