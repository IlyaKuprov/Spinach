# examples/microfluidics/reacting_flow_nmr.m

- Signature: `reacting_flow_nmr()`

## Purpose

Complete microfluidic simulation: diffusion, flow, two second- order chemical reactions, and NMR detection in a narrow strip of the chip where the coil is assumed to be located. Calculation time: days, much faster on GPU.

## Physical / mathematical content

- Microfluidics examples. The coupled model is spin dynamics plus advection-diffusion-reaction transport on a mesh or regular grid. Numerical issues include finite-difference operators, mesh interpolation, and coupled reaction-flow evolution.
- Signal processing is central here: the code moves between time and frequency domains, typically using FFT conventions, apodisation, zero filling, or heterodyne frequency shifts.
- Orientation or trajectory averaging is performed numerically, so grid design, weights, and integration error control matter directly to accuracy and runtime.

## Numerical / algorithmic content

- Time propagation is explicit. In Spinach this usually means repeated application of matrix exponentials or propagator factorizations to density operators or state vectors in Hilbert/Liouville/Fokker-Planck space.
- The output is processed in the Fourier domain, implying standard NMR/ESR signal-processing considerations such as acquisition bandwidth, zero filling, phase, and apodisation.
- The implementation explicitly addresses performance engineering through parallel or GPU execution, which matters because Spinach operators can become extremely large after basis expansion or powder/spatial lifting.
- Numerical integration over angles or geometry is part of the implementation, so point placement and weights are as important as the local Hamiltonian calculations.

## Implementation structure

- Complete microfluidic simulation: diffusion, flow, two second-
- order chemical reactions, and NMR detection in a narrow strip
- of the chip where the coil is assumed to be located.
- Calculation time: days, much faster on GPU.
- Import Diels-Alder cycloaddition
- Import hydrodynamics information
- Magnet field
- This needs a GPU
- Spinach housekeeping
- % Concentration dynamics stage
- Rate constants, mol/(L*s)
- Cycloaddition reaction generator, including solvent
