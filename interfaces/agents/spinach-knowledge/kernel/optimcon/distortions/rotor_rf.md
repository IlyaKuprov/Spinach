# Rotor-dependent RF gain

Bind a complex channel-by-tick array with @(w)rotor_rf(w,gain). Rows of w alternate X and Y. The map returns a sparse real Jacobian on column-major waveform samples. For one map per physical drift use dist_drift in ens_corrs; put electronics distortions before the sampled map.
