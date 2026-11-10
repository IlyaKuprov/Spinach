# rd_coupled_pair: radiation damping

[MATLAB example](../../../../../../examples/nmr_liquids/radiation_damping/rd_coupled_pair.m)

Compares quantum radiation damping of a weak AX doublet and a strongly coupled
AB pair with a classical model of their individual spectral lines. Both models
use one shared feedback field rather than a separate coil for each line.
The quantum calculation propagates the full spin density state; the line model
uses transition offsets and equilibrium intensities.

Longitudinal and transverse magnetisation traces expose the limitations of
replacing a coupled spin system by independent lines under nonlinear feedback.
The line representation is a comparison model, not a general substitute for
coupled-spin dynamics.
