# examples/kinetics/nonlinear/diels_alder_spec.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/examples/kinetics/nonlinear/diels_alder_spec.m`
- Signature: `diels_alder_spec()`
- Total lines: 265

## Purpose

Repeated pulse-acquire experiment during the Diels-Alder cyclo- addition of acetylene to butadiene, demonstrating the non-linear kinetics module. Calculation time: hours, GPU is hard-coded.

## Physical / mathematical content

- Chemical-kinetics examples. The files couple spin dynamics to exchange, pumping, or nonlinear reaction networks represented by kinetic generators in Liouville space.
- The relaxation model is Redfield-type perturbation theory: fluctuating interactions enter through correlation functions or spectral densities and generate a linear relaxation superoperator.
- Signal processing is central here: the code moves between time and frequency domains, typically using FFT conventions, apodisation, zero filling, or heterodyne frequency shifts.
- Orientation or trajectory averaging is performed numerically, so grid design, weights, and integration error control matter directly to accuracy and runtime.

## Numerical / algorithmic content

- Time propagation is explicit. In Spinach this usually means repeated application of matrix exponentials or propagator factorizations to density operators or state vectors in Hilbert/Liouville/Fokker-Planck space.
- The output is processed in the Fourier domain, implying standard NMR/ESR signal-processing considerations such as acquisition bandwidth, zero filling, phase, and apodisation.
- The implementation explicitly addresses performance engineering through parallel or GPU execution, which matters because Spinach operators can become extremely large after basis expansion or powder/spatial lifting.
- Numerical integration over angles or geometry is part of the implementation, so point placement and weights are as important as the local Hamiltonian calculations.

## Implementation structure

- Repeated pulse-acquire experiment during the Diels-Alder cyclo-
- addition of acetylene to butadiene, demonstrating the non-linear
- kinetics module.
- Calculation time: hours, GPU is hard-coded.
- DFT import options
- Load and display acetylene (substance A)
- Load and display butadiene (substance B)
- Load and display cyclohexadiene (substance C)
- Add natural abundance ethanol (substance D)
- Merge the spin systems
- Magnet field
- Greedy parallelisation

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `gparse()`, `g2spinach()`, `kfigure()`, `scale_figure()`, `subplot()`, `cst_display()`, `camorbit()`, `ktitle()`, `num2cell()`, `merge_inp()`, `create()`, `basis()`, `conc_traj()`, `step()`, `kxlabel()`, `kylabel()`.
