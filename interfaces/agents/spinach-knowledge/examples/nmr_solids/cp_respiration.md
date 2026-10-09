# examples/nmr_solids/cp_respiration.m

https://github.com/IlyaKuprov/Spinach/blob/main/examples/nmr_solids/cp_respiration.m

## Purpose

A 1H–13C RESPIRATION-CP example in the doubly rotating frame using magic-angle spinning and the Fokker–Planck formalism. The source estimates seconds of calculation time.

## Spin model and sequence inputs

The source sets sys.magnet to 11.7 and places 1H and 13C at [0, 0, 0] and [0, 0, 2.00]; these coordinates provide the pair geometry for the dipolar interaction. No additional interaction tensor is assigned in the wrapper. It uses the full sphten-liouv basis. The wrapper sets rate=20000 (no unit is stated), the axis [sqrt(2/3), 0, sqrt(1/3)], max_rank=8, and the rep_2ang_100pts_sph grid. It initialises 1H transverse magnetisation (Lx) and detects 13C L+; nloops is 16 and theta is pi/20. These are parameters passed to the external respiration sequence helper: this file does not define its internal RF waveform or give an RF-amplitude or CP-contact-time schedule.

## Acquired and plotted spectrum

The wrapper requests 512 points, zerofill to 16384, sweep=40000, and axis_units=kHz. It applies exponential apodisation with parameter 5, Fourier-transforms the FID, and plots the real spectrum. The code specifies the plotted frequency axis as kHz but does not annotate units for rate, sweep, or the apodisation parameter. This is a simulated output path, not a claim about a measured spectrum or a validation run.
