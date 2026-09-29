# kernel/grids/repulsion.m

[MATLAB source](https://github.com/IlyaKuprov/Spinach/blob/main/kernel/grids/repulsion.m) · [Spinach Wiki](https://spindynamics.org/wiki/index.php?title=repulsion.m) · [Bak–Nielsen reference](http://dx.doi.org/10.1006/jmre.1996.1087)

- Signature: `[alphas,betas,gammas,weights]=repulsion(npoints,ndims,niter)`
- Inputs: positive integer `npoints` and `niter`; `ndims` must be `2`, `3`, or `4`.
- Outputs: four `npoints x 1` columns. Angles are in radians; each `weights` entry is `1/npoints`, so the weights sum to one. They are uniform, not SHREWD-optimised.

The function starts with `rand(ndims,npoints)-0.5`, then repeats exactly `niter` updates. For each point pair it forms a normalised difference vector, sets non-finite entries from self-interactions to zero, and multiplies by the corresponding entries of `R'*R`; summing over partners gives `F`. It updates `R_new=R-ndims*F/npoints` and normalises every column to unit length. The initial random coordinates are not projected before the first update. The source prints the maximum displacement each iteration and does not stop on convergence.

Angle mapping depends on `ndims`:

- `2`: `betas=atan2(R(2,:),R(1,:))'`; `alphas` and `gammas` are zero columns. This is a one-angle circle grid.
- `3`: from `cart2sph`, `betas=theta+pi/2`, `gammas=phi`, and `alphas` is zero. With no output arguments, this branch also plots the grid.
- `4`: rows of `R` become the `u,i,j,k` components of quaternion records, which `qter2euler` converts to Euler angles.

The initial points depend on MATLAB's random-number state; the function does not seed the generator, so a fixed seed is needed for repeatable starts. It validates the stated scalar ranges but does not impose a convergence test or optimise the returned uniform weights. This function generates orientation grids; it does not calculate eigenfields, time evolution, or frequency offsets.
