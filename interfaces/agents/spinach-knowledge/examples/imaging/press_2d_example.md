# examples/imaging/press_2d_example.m

- MATLAB implementation: [examples/imaging/press_2d_example.m](https://github.com/IlyaKuprov/Spinach/blob/main/examples/imaging/press_2d_example.m)

## Purpose

This zero-argument example places three proton-pair components at separate positions in a two-dimensional phantom. It displays the phantom and the PRESS active-volume diagnostic, then runs a PRESS acquisition and plots its processed signal. The source estimates minutes of calculation time and notes that a Tesla V100 GPU makes it faster.

## Spin and spatial model

At 3.0 T, six 1H spins form three scalar-coupled pairs. Their assigned Zeeman scalar values are {-3, +3}, {-2, +2}, and {-1, +1}; the pair coupling values are 10, 20, and 30. The example does not state units for these values. It uses the sphten-liouv formalism, IK-2 approximation, proximity level 1, and scalar-coupling connectivity.

The sample dimensions are [0.30 0.25] on a [108 90] grid. Three circular phantom masks have radius 10 grid points and centres (15,15), (30,40), and (70,80). The image grid is [101 105]. Relaxation phantom values are zero; each phantom component has its own paired-spin Lz state. The receive-coil phantom is uniform and detects total L+.

The two slice-selection gradient amplitudes are 25e-3; the spatial-selection gradient amplitude is 5e-3 T/m and its duration is 5e-4 s. The PRESS callback documentation specifies T/m for slice-selection gradients, Hz for RF frequencies, rad/s for RF amplitudes, and seconds for pulse durations. The source does not label units for the sample dimensions.

## PRESS callback and readout

The source first plots the combined phantom, calls imaging with press_voxel_2d for an active-volume diagnostic, and then calls imaging with press_2d. The callback requires the imaging context to provide the evolution operators. The configured RF phases are {pi/2, pi/2}; the active frequency lists are {-120e3, -100e3}, amplitudes are {2*pi*5000, 2*pi*5000}, durations are {0.5e-4, 1.0e-4} s, and maximum ranks are {3, 3}.

The source comments list alternative RF frequencies for components B (`{-80e3, -10e3}`) and C (`{+30e3, +100e3}`), in Hz; only the A list `{-120e3, -100e3}` is active in this invocation. It applies square-cosine apodisation to the acquired data, takes a magnitude shifted Fourier transform, and plots the result with plot_1d. No numeric acquisition or fitted result is stated by the example.

## Source

[Example source](https://github.com/IlyaKuprov/Spinach/blob/main/examples/imaging/press_2d_example.m) · [PRESS callback](https://github.com/IlyaKuprov/Spinach/blob/main/experiments/imaging/press_2d.m) · [Voxel diagnostic](https://github.com/IlyaKuprov/Spinach/blob/main/experiments/imaging/press_voxel_2d.m)

The source credits Ahmed Allami and Ilya Kuprov.

Zero track elimination is explicitly enabled with `zte` in `sys.enable`.
