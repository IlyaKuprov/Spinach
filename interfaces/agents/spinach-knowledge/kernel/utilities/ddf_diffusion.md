# Reflecting voxel diffusion

`ddf_diffusion(mask,spacing,diffusivity)` couples only pairs of face-sharing interior voxels with the corresponding `D/h^2` rate. Columns sum to zero, conserving mass even in disconnected cavities; exterior voxels have no open faces. The operator is a rate, not a Hamiltonian. In the `exp(-1i*L*t)` convention its generator contribution is `1i*rate` tensored with the local spin identity. Voxel walls form a staircase approximation to a curved boundary; refine the mask independently of vacuum padding.
