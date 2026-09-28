# examples/microfluidics/reacting_nmr.m

- Signature: `reacting_nmr()`

## Purpose

Non-linear reaction kinetics in combination with spin evolution (repeated pulse-acquire NMR) and relaxation (Redfield theory). Calculation time: hours, much faster on GPU.

## Physical / mathematical content

- This is a homogeneous reaction-NMR calculation, not a spatial microfluidics model. It first integrates the two competing cycloaddition concentration kinetics, then couples the interpolated concentrations to chemical reaction superoperators during repeated proton pulse-acquire acquisitions.
- The relaxation model is Redfield-type perturbation theory: fluctuating interactions enter through correlation functions or spectral densities and generate a linear relaxation superoperator.
- Signal processing is central here: the code moves between time and frequency domains, typically using FFT conventions, apodisation, zero filling, or heterodyne frequency shifts.
- The script samples acquisitions at integer-second reaction times, constructs a waterfall spectrum, and uses two-point Lie quadrature for both the chemistry-coupled evolution and the NMR time steps.

## Numerical / algorithmic content

- The two reaction channels are integrated on a 20-second, 200-step grid. During each acquisition, the reaction-dependent generator is evaluated at both interval edges and passed to `step` with two-point Lie quadrature; relaxation is included in the spin generator.
- The output is processed in the Fourier domain, implying standard NMR/ESR signal-processing considerations such as acquisition bandwidth, zero filling, phase, and apodisation.
- The implementation explicitly addresses performance engineering through parallel or GPU execution, which matters because Spinach operators can become extremely large after basis expansion or powder/spatial lifting.
- The source enables greedy parallelisation and contains an optional GPU branch; it also apodises and zero-fills the assembled FIDs before the Fourier transform.

## Implementation structure

- Non-linear reaction kinetics in combination with spin evolution
- (repeated pulse-acquire NMR) and relaxation (Redfield theory).
- Calculation time: hours, much faster on GPU.
- Import Diels-Alder cycloaddition
- Magnet field
- Greedy parallelisation
- Spinach housekeeping
- Rate constants, mol/(L*s)
- Cycloaddition reaction generator, including solvent
- Kinetic time grid, 20 seconds
- Preallocate concentration trajectory
- Initial concentrations, mol/L
