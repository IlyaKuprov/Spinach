# examples/optimal_control/distortions/kernel_estimation/kernel_from_87rb.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/examples/optimal_control/distortions/kernel_estimation/kernel_from_87rb.m`
- Signature: `kernel_from_87rb()`
- Total lines: 130

## Purpose

Transmitter and probe distortion kernel of a 400 MHz Bruker spectrometer fitted with a 4 mm Phoenix MAS probe, estimated from an oscilloscope recording of an eight-block XiX wave- form on the 87Rb channel. The wall clock record is heterodyned into the rotating frame, the carrier frequency is refined from the residual phase drift accumulated within each XiX half-block, the remaining constant phase is removed, and the 

## Physical / mathematical content

- Optimal-control examples. These scripts formulate pulse design as a nonlinear optimisation problem over waveform samples or basis coefficients. The core mathematical objects are fidelities, gradients, Hessians or Hessian approximations, ensemble robustness objectives, and constrained search over RF amplitude/phase trajectories.
- Orientation or trajectory averaging is performed numerically, so grid design, weights, and integration error control matter directly to accuracy and runtime.

## Numerical / algorithmic content

- The output is processed in the Fourier domain, implying standard NMR/ESR signal-processing considerations such as acquisition bandwidth, zero filling, phase, and apodisation.
- Numerical integration over angles or geometry is part of the implementation, so point placement and weights are as important as the local Hamiltonian calculations.

## Implementation structure

- Transmitter and probe distortion kernel of a 400 MHz Bruker
- spectrometer fitted with a 4 mm Phoenix MAS probe, estimated
- from an oscilloscope recording of an eight-block XiX wave-
- form on the 87Rb channel.
- The wall clock record is heterodyned into the rotating frame,
- the carrier frequency is refined from the residual phase drift
- accumulated within each XiX half-block, the remaining constant
- phase is removed, and the FIR kernel is then obtained by line-
- ar least squares from the ideal and the measured waveforms.
- The heterodyne is an analytic signal demodulation, which is
- zero-phase: the rotating frame signal is not delayed with re-
- spect to the oscilloscope record.

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `load()`, `time_scope()`, `kfigure()`, `scale_figure()`, `subplot()`, `kxlabel()`, `kylabel()`, `ktitle()`, `heterodyne()`, `real_bruk()`, `imag_bruk()`, `n_low()`, `n_up()`, `polyfit()`, `time_grid()`, `unwrap()`.
