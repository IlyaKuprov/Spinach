# examples/nmr_overtone/dor_glycine.m

- Signature: `dor_glycine()`

## Purpose

Simulates a panoramic double-rotation (DOR) `14N` overtone spectrum of glycine, as described in Figure 1B of [the cited paper](https://doi.org/10.1039/C5CP03266K). The source describes a short pulse at instrumentally inaccessible power to make the excitation pattern uniform. The glycine quadrupolar tensor data are attributed to O'Dell and Ratcliffe, [Chemical Physics Letters (2011)](https://doi.org/10.1016/j.cplett.2011.08.030). The source estimates hours of calculation time.

## Model and calculation

The single-`14N` system is at 14.1 T, with quadrupole parameters 1.18 MHz and asymmetry 0.53 (spin 1). It uses diagonal damping relaxation at rate 100, the `sphten-liouv` basis without approximation, and disables Krylov and trajectory-level options.

At the magic angle, the experiment uses outer and inner spinning rates of 1425 and 6950 Hz, respectively, and rank 5 for each rotor. The outer axis is [sin(θ), 0, cos(θ)]; the inner axis is [√(20−2√30), 0, √(15+2√30)] (30.56°). The grid is `rep_2ang_200pts_sph`, the sweep is −20 to 30 kHz, and the acquisition uses 1024 points with 1024-point zero-fill. Average treatment uses RF power 2π × 3.0 MHz, a 1 μs pulse, and −10 kHz RF frequency. The code runs `doublerot` with `overtone_pa` and applies `exp(1i*1.49)` phase correction.
