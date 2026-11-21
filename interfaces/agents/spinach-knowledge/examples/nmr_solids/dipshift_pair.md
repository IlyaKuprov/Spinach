# examples/nmr_solids/dipshift_pair.m

`dipshift_pair()` compares signed DIPSHIFT echoes for a proton–carbon pair spinning at 6 kHz. At a fixed signed dipolar coupling of −10 kHz, the finite X,Y, zero-width X,Y, and X-only curves distinguish finite RF effects from changes in pulse schedule. The finite pulses are 5 μs pi pulses (100 kHz nutation frequency). An additional −23 kHz X-only curve shows the effect of changing coupling strength, not pulse width.

The example uses `singlerot` with Fokker–Planck MAS and a two-angle Lebedev powder grid. It plots the real echo divided by the prepared carbon overlap against pulse displacement in rotor periods. This scale retains the echo sign and does not force different curves to agree at zero displacement.
