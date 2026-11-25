# Tied Cartesian GRAPE blocks

`grape_repeat` repeats waveform columns in a full-duration `optimcon`
problem. Drift slices and dwell times remain independent at each time,
so this is not propagation by a repeated unitary. Gradients sum over
tied copies and Hessians contract both axes by the replication map.
Penalties and trajectories describe the expanded waveform. Waveform bases
and independent freeze masks cannot be used with this wrapper.
