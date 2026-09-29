# kernel/conventions/transforms/axis_tsymm.m

[Source](https://github.com/IlyaKuprov/Spinach/blob/main/kernel/conventions/transforms/axis_tsymm.m) · [Spinach Wiki](https://spindynamics.org/wiki/index.php?title=axis_tsymm.m)

## Contract

axis_tsymm forms an approximate rotational average of a real 3-by-3 interaction tensor about a specified Cartesian axis. It returns a 3-by-3 tensor in the same matrix representation; the routine does not take a spin system, orientation grid, or spatial-coordinate array.

The input axis is a nonzero real 3-by-1 vector, and the tensor is a real 3-by-3 matrix. The source does not specify a physical unit for the tensor entries or require the input matrix to be symmetric. The vector and tensor are interpreted in the common Cartesian frame used by the rotation operation.

## Transformation and sampling

For n=1:360, the routine constructs R=anax2dcm(a,pi*n/180) and accumulates R*T*R'; it divides the sum by 360. Thus the source samples 360 equally spaced rotations over a full turn, at pi*n/180 radians, and returns their discrete mean. “Approximate” describes this finite angular sampling.

## Source-supported use

The documented call is A=axis_tsymm(T,a), with T a 3-by-3 real interaction tensor and a the 3-by-1 real rotation-axis vector.
