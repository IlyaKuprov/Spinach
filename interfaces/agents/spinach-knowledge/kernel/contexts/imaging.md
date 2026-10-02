# kernel/contexts/imaging.m

- Signature: `answer=imaging(spin_system,pulse_sequence,parameters)`
- Source: [kernel/contexts/imaging.m](https://github.com/IlyaKuprov/Spinach/blob/main/kernel/contexts/imaging.m)
- Wiki: [Imaging.m](https://spindynamics.org/wiki/index.php?title=Imaging.m)

## Contract

`imaging` builds spin Hamiltonian and kinetic operators, voxel-dependent relaxation, spatial gradient operators, and the spatial diffusion/flow generator. It passes `H`, `R`, `K`, `G`, and `F` to `pulse_sequence(spin_system,parameters,H,R,K,G,F)`; the context output is whatever that sequence returns. The context applies the `nmr` assumption and channel-frequency offsets.

The spin-space dimension is `spn_dim=size(H,1)`; the spatial dimension is `spc_dim=prod(parameters.npts)`; the combined state-space dimension is their product. Spatial arrays are ordered as [X Y Z], while the direct-product factor order is Z, then Y, then X, then Spin. The corresponding vector is the column-wise vectorisation of a 3D [X Y Z] array with a spin-state component at each voxel. `parameters.spc_dim` and `parameters.spn_dim` are passed to the sequence.

## Grid, transport, and units

`parameters.dims` gives box lengths in metres and `parameters.npts` gives point counts along the axes. The flow fields `u`, `v`, and `w` are X-, Y-, and Z-velocity components at each sample point, in m/s. For spatially uniform diffusion, `parameters.diff` is a diffusion coefficient or 3-by-3 tensor in m^2/s. For voxel-dependent diffusion, the source documents Cartesian tensor-component fields from `dxx` through `dzz`, specified at each voxel.

The documented derivative choices are `parameters.deriv={'fourier'}` for Fourier differentiation, or `parameters.deriv={'period',n}` for n-point central finite differences with periodic boundary conditions.

With Fourier derivatives and `sys.enable={'polyadic'}`, the spatial generator `F` contains implicit FFT actions. The pulse sequence must use exponential-action propagation; `full(F)` and `inflate(F)` cannot materialise those cores. Disable polyadics before calling `imaging` when the sequence requires explicit spatial matrices or propagators. Finite-difference polyadics retain materialisable numeric cores.

## Spin operators and phantoms

The Hamiltonian and kinetics are shared across voxels. Relaxation is assembled from paired cell arrays `rlx_ph` and `rlx_op`: each `rlx_ph` entry is a spatial phantom with the same dimensions as the voxel grid, and the matching `rlx_op` entry is a spin relaxation superoperator. Initial states are assembled from `rho0_ph` and `rho0_st`; receiver/detection states are assembled from `coil_ph` and `coil_st`. Each phantom has the grid dimensions, and each paired state/operator is a spin-space object. The source requires the phantom and paired-state/operator lists to have matching lengths. A user-supplied `rho0` or `coil` can be used instead of building that object from phantoms.

## Example from the source documentation

For example, choose `parameters.deriv={'fourier'}` or `parameters.deriv={'period',n}` according to the requested spatial derivative scheme; the source gives these selector forms but no complete numerical imaging setup.
