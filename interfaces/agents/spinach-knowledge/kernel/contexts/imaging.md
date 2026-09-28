# kernel/contexts/imaging.m

- Signature: `answer=imaging(spin_system,pulse_sequence,parameters)`

## Purpose

Fokker-Planck imaging simulation context. Builds the Hamiltonian, relaxation and kinetics superoperators, spatial dynamics generator (including diffusion and flow), and gradient operators, then passes them to a pulse sequence supplied as a function handle.

## Physical / mathematical content

The context combines spin dynamics with a voxelized spatial state space. It applies the same Hamiltonian and kinetics superoperator in every voxel, while relaxation superoperators are weighted by voxel phantoms. The gradient and diffusion/flow generators are constructed by `g2fplanck` and `v2fplanck`, respectively. The direct-product order is `Z(x)Y(x)X(x)Spin`, corresponding to column-wise vectorization of a 3D array with dimensions `[X Y Z]`.

## Numerical / algorithmic content

The function checks parameter consistency, applies NMR assumptions, builds the Hamiltonian and processes channel offsets, then constructs the kinetics, relaxation, gradient, and spatial-dynamics operators. It builds initial and detection states from phantoms when explicit states are not supplied. Operators are represented as polyadic objects and inflated if `polyadic` is not enabled. Finally, it calls `pulse_sequence(spin_system,parameters,H,R,K,G,F)`.

## Parameters / inputs

- `pulse_sequence`: pulse-sequence function handle; see the Spinach `experiments` directory for supplied sequences.
- `parameters.u`, `parameters.v`, `parameters.w`: X, Y, and Z components of the velocity vectors at each sample point, in m/s.
- `parameters.diff`: diffusion coefficient or 3×3 tensor, in m²/s, when the diffusion parameter is the same in every voxel.
- `parameters.dxx`, `parameters.dxy`, …, `parameters.dzz`: Cartesian components of the diffusion tensor for each voxel.
- `parameters.dims`: dimensions of the 3D box, in meters.
- `parameters.npts`: number of points in each dimension of the 3D box.
- `parameters.deriv`: `{'fourier'}` selects Fourier differentiation matrices; `{'period',n}` selects n-point central finite-difference matrices with periodic boundary conditions.
- `parameters.rlx_ph={Ph1,Ph2,...,PhN}` and `parameters.rlx_op={R1,R2,...,RN}`: **required** relaxation phantom coefficients and corresponding relaxation superoperators. Both must be cell arrays with the same number of elements; each phantom must match the voxel grid specified by `parameters.npts`.
- `parameters.rho0_ph={Ph1,Ph2,...,PhN}` and `parameters.rho0_st={rho1,rho2,...,rhoN}`: initial-condition phantoms and corresponding spin states from `state()`. These are required **only when `parameters.rho0` is absent**; otherwise the supplied initial state is used.
- `parameters.coil_ph={Ph1,Ph2,...,PhN}` and `parameters.coil_st={rho1,rho2,...,rhoN}`: detection-state phantoms and corresponding spin states from `state()`, allowing voxel-dependent detection angles and sensitivities. These are required **only when `parameters.coil` is absent**; otherwise the supplied coil state is used.

Each initial-condition or detection phantom used must match `parameters.npts` and the relaxation phantoms in size. The corresponding phantom and state lists must be cell arrays of equal length.

## Outputs

Returns whatever the pulse sequence returns.

## Reference

- <https://spindynamics.org/wiki/index.php?title=Imaging.m>

## Authors

- a.j.allami@soton.ac.uk
- ilya.kuprov@weizmann.ac.il