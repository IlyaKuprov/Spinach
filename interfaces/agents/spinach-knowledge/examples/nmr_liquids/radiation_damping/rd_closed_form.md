# rd_closed_form: radiation damping

[MATLAB example](../../../../../../examples/nmr_liquids/radiation_damping/rd_closed_form.m)

Integrates radiation damping of a proton initially tilted almost opposite to
its equilibrium magnetisation. The state-dependent feedback generator is
propagated with fourth-order Lie-group product quadrature, without relaxation.
The longitudinal recovery and transverse envelope are compared directly with
the Bloembergen–Pound tanh and sech solutions.

Run this example first to check the feedback sign and equilibrium normalisation.
The closed form applies to a single resonant spin, not to a dephasing ensemble
or a coupled multiplet.
