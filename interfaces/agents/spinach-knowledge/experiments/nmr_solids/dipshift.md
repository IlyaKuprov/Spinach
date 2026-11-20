# experiments/nmr_solids/dipshift.m

`signal=dipshift(spin_system,parameters,H,R,K)` is a Liouville-space callback for `singlerot`. Set `parameters.spins={recoupled,observed}` to two distinct isotope labels, `rate` to the signed MAS rate in Hz, and `time_shift` to a nonempty row of displacements in seconds from zero to half a rotor period. Supply the prepared state and receiver in `rho0` and `coil`; observed-channel `Lx` preparation and `L+` detection give a negative echo without coupling.

Set `xy_block=true` for one X,Y recoupling block in each arm, followed or preceded by the central echo. The echo drives X on the recoupled channel and Y on the observed channel simultaneously. Set `pulse_width` in seconds: positive widths include drift during the pulse, and zero width keeps the same pulse centres and pi rotations. The alternative `xy_block=false` is an X-only schedule with different timing and requires zero pulse width. Do not use it as the hard-pulse limit of the X,Y sequence.

The result is an unnormalised complex row in displacement order, not a directly sampled FID. Each point starts from the same prepared state and rotor distribution. Rotor-rank and powder-grid convergence are separate from pulse-width convergence.
