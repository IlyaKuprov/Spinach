# examples/optimal_control/distortions/kernel_estimation/kernel_from_87rb.m

- MATLAB implementation: [examples/optimal_control/distortions/kernel_estimation/kernel_from_87rb.m](https://github.com/IlyaKuprov/Spinach/blob/main/examples/optimal_control/distortions/kernel_estimation/kernel_from_87rb.m)

- Signature: kernel_from_87rb()
- Source: [examples/optimal_control/distortions/kernel_estimation/kernel_from_87rb.m](../../../../../../../examples/optimal_control/distortions/kernel_estimation/kernel_from_87rb.m)

## Purpose and recording

Estimates a complex transmitter/probe finite-impulse-response distortion kernel from an oscilloscope recording of an eight-block XiX waveform on the 87Rb channel. The source describes the setup as a 400 MHz Bruker spectrometer with a 4 mm Phoenix MAS probe. It loads 87Rb_160W_XiX_10us_8repeats.mat (scope_trace, time_scope); the recorded trace is imported data, not a Spinach simulation.

## Carrier and phase processing

The requested kernel grid is 100 ns with 40 taps (4 μs at that sampling interval). The XiX train has eight 10 μs blocks, starts 11.1 μs into the oscilloscope record, and has a nominal 87Rb carrier of 131.1468520 MHz. Fitting windows occupy fractions 0.32–0.48 of each 5 μs half-block, away from transitions. The code heterodynes the record at the nominal carrier, unwraps the phase within each window, and uses a linear polyfit slope for each of the 16 half-blocks. It refines the carrier by subtracting the mean residual slope from all but the first window, then heterodynes again. A constant phase is removed using the mean phase of the nominally positive half-block windows, and the complex segment is peak-normalised before resampling onto the 100 ns grid.

## FIR estimate and outputs

The reference is an ideal square XiX envelope that alternates sign every half-block and is zero after the 80 μs train. kernelest fits a 40-tap complex kernel to this reference and the processed oscilloscope record using the pinv method; dividing by sum(h) normalises its DC component. The figure shows the raw wall-clock trace, heterodyned in-phase/quadrature signal against the ideal waveform, complex impulse-response components over the 40-tap time support, and the zero-filled magnitude frequency response (256-point FFT). The script saves h and dt_kernel in kernel_87rb.mat.

This documents the source's signal-processing procedure and reported instrument context, not an independently checked acquisition or validation. The actual oscilloscope MAT file is required to reproduce the estimate; the script provides no uncertainty analysis or separate kernel-validation metric.
