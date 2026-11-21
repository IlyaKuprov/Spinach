# experiments/nmr_solids/respdor.m

Use `respdor` as a Liouville-space `singlerot` callback for an observed spin-1/2 isotope and a distinct quadrupolar isotope, in that order. The R schedule uses two rotor periods of rotary resonance per arm; the S schedule uses a six-period phase-supercycled R4 block. The second arm repeats the first in the same order.

Both echoes span two rotor periods between arms. The finite observed-channel inversion sits at the centre of the irradiation window. The quadrupolar RF is omitted only in the reference; drift and elapsed time are retained, including during the inversion. Zero cycles report the unrefocused prepared-state overlap.

The output rows contain complex reference, irradiated, and reference-minus-irradiated signals. Average these rows coherently before taking their ratio; omit points with unresolved reference signal.
