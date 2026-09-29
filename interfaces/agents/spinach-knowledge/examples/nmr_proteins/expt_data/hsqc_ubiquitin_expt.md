# examples/nmr_proteins/expt_data/hsqc_ubiquitin_expt.m

Source: [examples/nmr_proteins/expt_data/hsqc_ubiquitin_expt.m](https://github.com/IlyaKuprov/Spinach/blob/main/examples/nmr_proteins/expt_data/hsqc_ubiquitin_expt.m)

## Purpose

Processes and plots measured two-dimensional 15N-1H HSQC data for human ubiquitin. The function loads `fid` from `hsqc_ubiquitin_expt.mat`; it does not simulate a pulse sequence, and its `spin_system` struct is used only for plotting metadata. The source specifies no paramagnetic centres or magnetic tensors, so this is protein nuclear-spin NMR rather than a paramagnetic calculation.

## Data processing

The loaded FID has positive and negative components. Each component is cosine-apodised in both dimensions and multiplied by the F1 phase factor `exp(-1i*0.7)`. The code Fourier-transforms each component along dimension 1 with zero filling to 1024 points, combines them as `f1_pos + conj(f1_neg)` to form a States signal, then Fourier-transforms dimension 2 to 1024 points. It flips both dimensions and plots `-imag(spectrum)`; it does not write a processed spectrum file.

## Axes and display

Plotting metadata specifies spins `15N`, `1H`, field `11.7395` T, offsets `[-5870 3753]` Hz, sweeps `[2000 4000]` Hz, and zero filling `[1024 1024]`; axes are displayed in ppm. The source processes imported experimental data and does not specify the pulse sequence or its receiver acquisition.
