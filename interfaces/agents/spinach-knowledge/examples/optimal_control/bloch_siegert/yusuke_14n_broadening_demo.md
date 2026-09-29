# examples/optimal_control/bloch_siegert/yusuke_14n_broadening_demo.m

> **Historical example — no longer shipped in current Spinach.** This page describes a deleted example only; do not try to run it or assume its source path exists in the current package.

- Historical function signature: yusuke_14n_broadening_demo()

## What the historical script illustrated

The script was a qualitative reduced-model illustration of how 14N decoupling conditions can affect an observed 1H line. Its comment describes a magnetic field corresponding to 800 MHz 1H, and the code sets magnet=18.8. It contrasted low-power CW, high-power CW, and a manually specified low-power offset-tolerant profile. That last profile was an illustrative design target, not a waveform produced by optimal control. The model was not a full quadrupolar/MAS calculation or a quantitative prediction.

It used three 14N site offsets (−18, 0, +18 kHz), weighted 0.30/0.40/0.30, and an 81-point B1 scale grid from 0.85 to 1.15 with Gaussian weighting arguments 1.0 and 0.06. The scenario RF values were 8 kHz (low-power CW), 20 kHz (high-power CW), and 12 kHz (offset-tolerant profile). The line-shape model also set an intrinsic-width parameter of 80, a residual-penalty parameter of 220, and a Bloch–Siegert coefficient of 8×10⁻⁷; these are model inputs, not reported experimental measurements. In the toy calculation, the offset response is eta = 1/(1 + (|offset|/(bandwidth_gain × RF))^profile_order), and each site line width is intrinsic width + residual penalty × (1 − eta). The Bloch–Siegert shift scales as coefficient × (RF × B1 scale)^2; after subtracting the ensemble-mean shift, the script sums weighted Lorentzian lines over the site and B1 distributions.

The script forms predicted 1H line shapes and calculates mean decoupling efficiency and numerical FWHM for the scenarios, including homogeneous and distributed B1 cases for high-power CW. Its intended qualitative point is that B1 dispersion converts a drive-induced Bloch–Siegert shift into line broadening. The source contains calculation and print/plot statements, but this note does not claim run-derived numerical outcomes.

The source names Nehra, Agarwal, and Nishiyama as inspiration and calls the manually specified profile a design target for Bloch–Siegert-aware robust optimal control, but supplies no complete reference or DOI; none is inferred here.
