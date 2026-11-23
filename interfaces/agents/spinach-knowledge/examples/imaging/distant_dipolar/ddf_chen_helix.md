# Transverse-helix frequency shift in spherical water

`ddf_chen_helix()` implements the ideal bulk limit of the matched-gradient
experiment of Chen, Branca, and Warren, [Fig. 1(b), Fig. 4 and Eq. (22)](https://doi.org/10.1063/1.4721637).
The spherical demagnetisation factor sets the uniform field to zero.
The weak transverse helix nonetheless generates a transverse field that
torques the large uniform longitudinal magnetisation. A longitudinal-only
source misses that frequency shift entirely.

The decoding kick is represented by the adjoint spatial phase in the
receiver, so one trajectory supplies all evolution delays. Phase slopes
come from the complex ratio to the dipolar-off echo, not from the magnitude
spectrum or the strength of the transverse field. Finite pulse durations,
RF refocusing, diffusion, relaxation, radiation damping, and finite sphere
edges are not modelled. The example therefore checks the ideal phase shift,
not the experimental envelope.

With the reported dipolar time of 130 ms the ideal axial phase slope has
magnitude `abs(cos(tip))/(3*tau_d)`, about 2.55 rad/s (0.41 Hz) for 5°.
This compares with the paper's rounded Eq. (22) estimate of 2.5 rad/s,
described there as close to the measured shift. Figure 4 shows the axial
shift and magic-angle null for 5°; the 175° reversal is reported in the
text. The 0.5° calculation probes the finite small-tip limit and is not an
additional experimental point. All reported signs use the raw Spinach
receiver convention.
