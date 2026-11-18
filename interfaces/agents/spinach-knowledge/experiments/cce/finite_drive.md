# finite_drive

Syntax: `[answer,aux]=finite_drive(spin_system,parameters,H,R,K)`.

A native five-argument finite-drive experiment for clexpn. Each element of the nonempty
parameters.drive structure row requires spins (observed physical IDs), axis (x/y/z or a
real unit three-vector), amplitude and detuning in rad/s, and positive duration
in seconds. Equal-length amplitude/duration rows specify piecewise-constant
waveform slices. Optional time is a scalar start (default zero), or a matching
row of explicit slice starts. Without explicit starts, slices are consecutive.
Observations may occur before, during, or after a pulse; propagation stops at each requested time. Alternatively, a single finite pulse may specify fraction in [0,1] instead of time; its duration remains in seconds. Multiple drive descriptors may have different axes, starts, and detunings, and overlapping amplitudes add. Drift remains active throughout.

Optional refocusing is a cell row of ideal-pulse descriptors: spins, either
an integrated Hermitian operator or axis and angle, and either time (seconds)
or fraction of the observation time. Coincident pulses keep declaration order.
This mixed ideal/finite path requires generalised Hilbert dynamics. Pure finite
drives retain the existing generalised master-equation path. A one-site exterior
pulse is inactive without mean fields; finite exterior excitation with static
mean fields and arbitrary exterior operator rotations are refused.

Preparation, detection, normalisation, density layout, and dissipative generators
follow free_evolution. Projected generators are unsupported. Retained bath drift
stays active throughout the pulse. The Hilbert-space finite-control path uses
segmentwise expm transport and its existing arithmetic bound; it does not invoke
the legacy driver. cce_acquire compiles native declarations into internal events. A caller-supplied sequence field is rejected even when empty.

Example: gcce_finite_drive. Validation: test_clexpn_experiments and the independent
dense reference in test_cce_examples.

Elapsed times and native drive/ideal-control descriptor fields are checked locally before acquisition; validation does not compile and discard a pulse-event array.
