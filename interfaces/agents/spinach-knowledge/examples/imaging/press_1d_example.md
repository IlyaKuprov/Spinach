# examples/imaging/press_1d_example.m

- MATLAB implementation: [examples/imaging/press_1d_example.m](https://github.com/IlyaKuprov/Spinach/blob/main/examples/imaging/press_1d_example.m)

## Purpose

This zero-argument example configures a one-dimensional PRESS localisation experiment for three spatially distinct proton-pair components, then displays a voxel-selection diagnostic and a spectrum. The example source estimates seconds of calculation time and notes that a Tesla V100 GPU makes it faster.

## Spin and spatial model

The model is at 3.0 T and contains six 1H spins arranged as three scalar-coupled pairs. The source assigns Zeeman scalar values {-3, +3}, {-2, +2}, and {-1, +1}, with pair coupling values 10, 20, and 30; it does not label units for those numbers. The basis is sphten-liouv with IK-2, proximity level 1, and scalar-coupling connectivity. Path tracing is disabled; GPU enablement is shown only as a commented option.

The sample extent is 0.30 and its one-dimensional phantom has 100 points. Three initial-state phantom profiles are supplied, with matching pair-spin states; the receive coil has a uniform spatial profile and detects total L+. The example sets both slice-selection and readout gradient amplitudes to 30e-3; the PRESS callback documentation specifies T/m for the slice-selection gradient. The source does not label the sample-extent unit here.

## PRESS callback and readout

The example first calls imaging with press_voxel_1d to display the excitation profile. It then calls imaging with press_1d, the imaging-context sequence callback. The PRESS callback documentation defines RF frequencies in Hz, RF amplitudes in rad/s, pulse durations in seconds, and requires the imaging context to supply the evolution operators. The configured RF phase is pi/2, frequency is -100e3, amplitude is 2*pi*5000, duration is 0.5e-4 s, and maximum rank is 3.

The source comments describe possible excitation assignments at +100e3 for B and C, 0 for C, and -100e3 for A and C, but the active example sets only -100e3; it does not execute a frequency scan. After acquisition, the example applies square-cosine apodisation, computes the magnitude of a shifted Fourier transform, and plots the result as a voxel spectrum. These are the requested processing steps, not reported numerical simulation results.

## Source

[Example source](https://github.com/IlyaKuprov/Spinach/blob/main/examples/imaging/press_1d_example.m) · [PRESS callback](https://github.com/IlyaKuprov/Spinach/blob/main/experiments/imaging/press_1d.m) · [Voxel diagnostic](https://github.com/IlyaKuprov/Spinach/blob/main/experiments/imaging/press_voxel_1d.m)

The source credits Ahmed Allami and Ilya Kuprov.
