# examples/nmr_solids/hmqc_mas_sq.m

- Signature: `hmqc_mas_sq()`

## Purpose

Simulates a rotor-synchronised powder MAS CN2D experiment on a `14N`–`1H` pair using a one-dimensional Fokker–Planck treatment and spherical orientation grid. The model includes the second-order quadrupolar shift and lineshape. The source estimates hours on CPU or minutes on a Tesla V100 GPU.

## Physical and numerical content

The system is specified at 19.96 T with a `14N` quadrupolar interaction (`Cq = 1.18 MHz`, `eta = 0.53`, spin `1`) and a `1H` spin. The MAS rate is 125 kHz about `[1 1 1]`; the two-dimensional sweep widths are 31.25 and 20 kHz. The powder grid is `rep_2ang_200pts_sph`. The experiment uses the `cn2d_sq` sequence. The cosine and sine FIDs are apodised, Fourier transformed in both dimensions, combined as a States signal, and plotted.

## Implementation

The function builds an exact sphten-liouv basis (maximum rank 16), calls `singlerot` for the SQ CN2D acquisition, then performs the two-dimensional processing. The source sets 256 × 128 acquired points and zero-fills to 1024 × 512; the `14N` rotating frame is specified at rank 3. RF power and duration are 40 kHz and 2 ms.
