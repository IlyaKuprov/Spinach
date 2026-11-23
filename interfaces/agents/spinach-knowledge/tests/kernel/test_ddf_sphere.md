# Sphere diffusion and bulk field checks

The radial Neumann mode `sin(q*r)/(q*r)` in a unit sphere, with `tan(q)=q`, decays to `exp(-1)` of its initial amplitude after `1/q^2` at unit diffusivity. Two voxel masks test convergence of the full radial mode while checking conservation and exterior leakage. The uniform-sphere cavity field is checked inside half the radius to separate bulk behaviour from surface voxelisation.
