# cce_segments

[MATLAB source](../../../../../kernel/cce/cce_segments.m)

Internal syntax: `segments=cce_segments(system,parameters,drift,events)`

Compiles experiment-owned events for one native retained subsystem. Physical
spin order, pulse support order, and simultaneous-event order are retained.
Intervals split at pulse starts/stops and observation endpoints, including the
initial zero-duration interval. Returned durations, Hamiltonians, pulse matrices,
and post-pulse exterior expectations expose the actual numerical calculation.

The segment engine reuses generator-supplied Cartesian spin operators. Finite
pulses add to the full drift. Static exterior Cartesian rotations update the
same bath product densities and subsequent fields through cce_subsystem and
cce_generators; finite exterior excitation or explicit exterior operator pulses
are refused. Integrated generators must be Hermitian relative to their own norm.
Retained segment arrays and propagation workspaces are counted before allocation.
No legacy model/parameter translation or propagation-result wrapper is provided.

Caller-supplied parameters.sequence is rejected, including empty arrays. Use finite_drive, refocusing_train, or the appropriate named echo experiment; internal events are constructed from native physical declarations, never read from that removed field.

Observable detection may omit coil when a transition supplies the detector. Segment memory accounting retains one detection channel in that case; ideal trains and mixed finite/ideal drives use the same optional-field contract.

Selective operator pulses are embedded directly in retained tensor order without changing their declared support ordering. Exterior-support restrictions are applied only to pulses actually included at an observation endpoint; omitted future pulses do not acquire new execution requirements.
