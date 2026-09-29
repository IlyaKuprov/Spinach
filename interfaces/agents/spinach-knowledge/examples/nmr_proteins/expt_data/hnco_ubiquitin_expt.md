# examples/nmr_proteins/expt_data/hnco_ubiquitin_expt.m

Source: [examples/nmr_proteins/expt_data/hnco_ubiquitin_expt.m](https://github.com/IlyaKuprov/Spinach/blob/main/examples/nmr_proteins/expt_data/hnco_ubiquitin_expt.m)

## Purpose

Processes and plots a three-dimensional experimental HNCO spectrum of human ubiquitin. This is processing of imported measured data, not a Spinach spin-dynamics simulation: the function loads `fid` from `hnco_ubiquitin_expt.mat`, does not call a simulator or pulse-sequence function, and uses its `spin_system` struct only for plotting metadata. The source does not specify paramagnetic centres or magnetic tensors; this is protein nuclear-spin NMR.

## Data processing

The code truncates the first FID dimension to 64 points and applies cosine apodisation in all three dimensions. It Fourier-transforms F3 to 256 points and applies a phase factor of `exp(-1i*0.75)`. For F2 it recombines alternating real components as `real(odd) + 1i*real(even)` before a 256-point transform; F1 uses `real(odd) - 1i*real(even)` before its 256-point transform. It permutes dimensions to `[3 2 1]`, shifts each dimension, subtracts `spectrum(end,end,end)` as a baseline, and zeros indices 1 through 40 in the third array dimension. The result is plotted as its real part; no output spectrum file is written.

## Axes and display

The plotting metadata specifies spins `15N`, `13C`, `1H`, magnetic field `11.7395` T, offsets `[-5870 22164 3653]` Hz, sweeps `[2000 1500 3000]` Hz, 64 acquired points per dimension, and zero filling to `[256 256 256]`. The axis display is in ppm. The source does not encode the experimental pulse sequence or receiver-channel acquisition; it processes the already acquired FID.
