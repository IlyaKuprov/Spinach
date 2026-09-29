# examples/nmr_solids/case_studies/cp_square_vs_ramp/cp_adiabatic_vs_optimcon.m

- MATLAB implementation: [examples/nmr_solids/case_studies/cp_square_vs_ramp/cp_adiabatic_vs_optimcon.m](https://github.com/IlyaKuprov/Spinach/blob/main/examples/nmr_solids/case_studies/cp_square_vs_ramp/cp_adiabatic_vs_optimcon.m)

Source: https://github.com/IlyaKuprov/Spinach/blob/main/examples/nmr_solids/case_studies/cp_square_vs_ramp/cp_adiabatic_vs_optimcon.m

This example compares simulated ¹H–¹⁵N cross-polarisation in the doubly rotating frame using a tangent-ramped contact and GRAPE-optimised controls. The spin pair is ¹⁵N and ¹H; the source sets the field parameter to 9.394 and temperature to 298, without annotating units for those values. The coordinates are [0 0 0] and [0 0 1.05], and both scalar Zeeman inputs are zero; coordinate units are not stated. No experimental dataset is loaded: the powder-averaged trajectories and waveforms are simulation products.

The common experiment uses 500 time steps of 2e-6 seconds each, for a 1 ms contact. The irradiation operators are Ly on ¹H and Lx on ¹⁵N; the excitation operator is Lx on ¹H and the detected operator is Lx on ¹⁵N. The powder grid is rep_2ang_200pts_sph. The tangent-ramped contact builds tan(linspace(-1.4, 1.4, 500)), normalises it to 0–1, and applies the reversed and forward ramps to the two channels, with a peak nutation-frequency input of 5e4 Hz. The source plots channel nutation frequencies in Hz against time in seconds.

For the GRAPE comparison, the initial state is normalised Ly on ¹H and the target is normalised Lx on ¹⁵N. The two control channels use a power level of 2π × 5e4 rad/s, the same 500 slice durations, an SNS penalty with weight 100, L-BFGS optimisation, and a 30-iteration limit. The initial guess is the tangent-ramp waveform. A second optimisation halves each slice duration, giving a 0.5 ms contact while retaining the same number of slices. The resulting controls are simulated with the same powder-averaged CP routine.

The figure compares the two-channel nutation-frequency waveforms and the real ¹⁵N X-expectation trajectory for the tangent ramp, the same-duration GRAPE pulse, and the half-duration GRAPE pulse. The code does not report an experimental transfer measurement or a quantitative agreement statistic.
