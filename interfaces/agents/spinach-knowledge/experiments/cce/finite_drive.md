# finite_drive

Syntax: `[answer,aux]=finite_drive(spin_system,parameters,H,R,K)`.

A native five-argument rectangular-drive experiment for clexpn. The scalar
parameters.drive descriptor contains spins (unique physical observed-spin IDs),
axis (x/y/z or a real unit three-vector), amplitude and detuning (real finite
rad/s scalars), and duration (positive finite seconds). All five fields are
required; unknown fields and a competing sequence array are rejected. The drive
starts at zero. Observation times are zero or at least the drive duration.

Preparation, detection, normalization, density layout and dissipative generators
follow free_evolution. Projected generators are unsupported. Retained bath drift
stays active throughout the pulse. The Hilbert-space finite-control path uses
segmentwise expm transport and its existing arithmetic bound; it does not invoke
the legacy driver. Pulse-event construction belongs to the experiment.

Example: gcce_finite_drive. Validation: test_clexpn_experiments and the independent
dense reference in test_cce_examples.
