# kernel/conventions/transforms/zfs2mat.m

Source: [kernel/conventions/transforms/zfs2mat.m](https://github.com/IlyaKuprov/Spinach/blob/main/kernel/conventions/transforms/zfs2mat.m)
Wiki: [Spin Dynamics Wiki: zfs2mat.m](https://spindynamics.org/wiki/index.php?title=zfs2mat.m)

- Signature: `M=zfs2mat(D,E,alp,bet,gam)`

## Purpose and tensor convention

Converts the zero-field-splitting parameters `D` and `E` into the symmetric spin-interaction matrix used by Spinach. In the tensor eigenframe the source first forms the diagonal matrix with entries `-D/3+E`, `-D/3-E`, and `2*D/3`. It computes the direction-cosine matrix `R=euler2dcm(alp,bet,gam)` and rotates the tensor as `M=R*M*R'`. It then removes any residual trace and symmetrises the result. The trace correction and symmetrisation are explicit numerical clean-up steps; the function does not solve an eigenproblem or compute a numerical derivative.

## Inputs and output

- `D`, `E` — real numeric scalar zero-field-splitting parameters in Hz.
- `alp`, `bet`, `gam` — real numeric scalar Euler angles in radians.
- `M` — symmetric `3x3` interaction matrix in Hz.

The implementation's `grumble` guard rejects any input that is not a real numeric scalar, with the message “all inputs must be real scalars.” There are no other input-dependent branches in this routine.

## Reference

The source cites the zero-field-splitting convention in the abstract of [doi:10.1063/1.1682294](http://dx.doi.org/10.1063/1.1682294). The source and Wiki describe the operation but provide no worked numerical example; none is added here.
