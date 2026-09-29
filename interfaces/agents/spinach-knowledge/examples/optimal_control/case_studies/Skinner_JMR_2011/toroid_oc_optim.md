# examples/optimal_control/case_studies/Skinner_JMR_2011/toroid_oc_optim.m

[Source](https://github.com/IlyaKuprov/Spinach/blob/main/examples/optimal_control/case_studies/Skinner_JMR_2011/toroid_oc_optim.m) · [Skinner et al., J. Magn. Reson. 209, 282–290 (2011)](https://doi.org/10.1016/j.jmr.2011.01.026)

Designs a *new Spinach* optimal-control excitation pulse under the paper's toroid geometry, not a copy of the article's original optimised waveform. The radial positions 1–6 mm are sampled uniformly so that the outer-wall 25 kHz RF amplitude scales as inverse radius; nine offsets cover ±1.5 kHz. The proton starts along z and is driven toward x over a 22 µs pulse digitised into 0.5 µs steps. The waveform is optimised by L-BFGS over a radius–offset ensemble with a soft amplitude spillout penalty, then rescaled if necessary to respect the 25 kHz outer-wall peak.

`[result,fig]=toroid_oc_optim()` returns the physical quadratures in Hz, the denser 51×51 verification response, and the paper's independently simulated 3.74 µs hard-pulse baseline. The plot shows the shaped-pulse offset–radius response and detected radius-weighted profile against that baseline. An optimisation result must be validated numerically before its performance is cited; its design is a source-geometry analogue rather than an exact recreation of a published pulse.
