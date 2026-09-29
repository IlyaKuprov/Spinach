# examples/optimal_control/case_studies/Kobzar_JMR_2004/bebop_bibop_profiles.m

[Source](https://github.com/IlyaKuprov/Spinach/blob/main/examples/optimal_control/case_studies/Kobzar_JMR_2004/bebop_bibop_profiles.m) · [Kobzar et al., J. Magn. Reson. 170, 236–243 (2004)](https://doi.org/10.1016/j.jmr.2004.06.017)

Simulates two original author-deposited single-spin Cartesian waveforms: a 337.5 µs BEBOP excitation pulse (longitudinal to transverse x magnetization) and a 307.5 µs BIBOP inversion pulse (longitudinal magnetization reversed). The [excitation](https://www.ioc.kit.edu/luy/186.php) and [inversion](https://www.ioc.kit.edu/luy/227.php) files contain RF quadratures in Hz and 0.5 µs slice durations in seconds. Their peak RF amplitudes are 10 kHz.

Calling `[profiles,fig]=bebop_bibop_profiles()` propagates both source waveforms in Spinach over the paper's 20 kHz total offset span (200 points) and five RF scalings from 0.8 to 1.2. The returned structure contains two offset-by-RF transfer-efficiency arrays, and the figure shows the original quadratures and independent two-dimensional response maps. Optional offset and RF-scale vectors permit smaller checks without changing the defaults. The arrays are state-transfer scores, not a claim to reproduce experimental data.
