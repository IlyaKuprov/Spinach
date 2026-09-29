# experiments/imaging/spiral.m

[MATLAB source](https://github.com/IlyaKuprov/Spinach/blob/main/experiments/imaging/spiral.m) · [Spinach Wiki](https://spindynamics.org/wiki/index.php?title=spiral.m)

## Purpose

`spiral` simulates a two-dimensional MRI acquisition with spiral k-space sampling. In the imaging context it applies 90° excitation, evolves for `t_echo`, applies a 180° pulse, evolves for the second `t_echo`, and then traverses a parameterised spiral while recording the selected observable. It constructs an image computationally; the source does not report a measured scan.

## Inputs and units

The function receives `H`, `R`, `K`, `G`, and `F` from `imaging()` and forms `L=H+F+1i*R+1i*K`. Its sequence-specific fields are:

- `parameters.t_echo`: echo time in seconds, used for each of the two free-evolution delays.
- `parameters.spiral_frq`: spiral angular frequency in rad/s.
- `parameters.spiral_dur`: spiral duration in seconds.
- `parameters.spiral_npts`: number of time samples along the trajectory.
- `parameters.grad_amp`: gradient amplitude at the end of the spiral, in T/m.
- `parameters.spins`: cell array of spin-name strings; the first entry selects the pulse operator and gyromagnetic ratio.
- `parameters.rho0`, `parameters.coil`, and `parameters.npts`: initial state, detection observable, and imaging-grid dimensions provided through the imaging setup.

The trajectory is formed from in-plane `G{1}` and `G{2}` operators. The code accumulates `spiral_kx` and `spiral_ky` in cycles per metre using the selected spin's gyromagnetic ratio; it records `parameters.coil'*rho` before each propagation step. The pulse and readout coherence selection is therefore determined by the selected spin operator and supplied coil observable, rather than an explicit coherence-order filter.

## Detection and return

The acquired samples are interpolated onto a square k-space grid with cubic interpolation, NaN cells are set to zero, square-sine apodisation is applied in both dimensions, and a shifted two-dimensional Fourier transform produces `mri`. The return is an image array, not separate `x`/`y` coordinate vectors; the displayed sampling trajectory is labelled in cycles per metre. The square grid uses `sqrt(spiral_npts)/2` points along each dimension, as coded in the source. No microfluidic-flow model or measured flow result is specified by this function.
