# rd_inhomogeneous: radiation damping

[MATLAB example](../../../../../../examples/nmr_liquids/radiation_damping/rd_inhomogeneous.m)

Propagates Gaussian-weighted proton isochromats in one shared detection coil.
Their weighted transverse magnetisation generates the same feedback field on
every isochromat; evolving independent feedback problems and averaging afterwards
would describe a different physical system. Plots compare longitudinal recovery
and transverse decay at several Gaussian full widths at half maximum, including
a resonant zero-width line.

The symmetric grid contains zero exactly. Its range reaches a relative Gaussian
height of 1e-6 for the widest line, and its spacing puts the discrete recurrence
at least four acquisition durations away. Recompute this grid when changing the
linewidths or acquisition duration, and check refinement for new parameters.
The calculation has no phenomenological relaxation.
