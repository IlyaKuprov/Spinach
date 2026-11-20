# examples/nmr_solids/respdor_pair.m

Run `respdor_pair()` to compare R- and S-RESPDOR on the same model 13C–14N pair at 10 kHz MAS, zero offsets, and 40 kHz quadrupolar RF. The schedules retain their distinct preparation, detection, recoupling, and irradiation timing. The horizontal axis adds the durations of the two recoupling arms, excluding the common echo block.

The example uses reduced quadrupolar coupling and a modest powder grid to demonstrate the sequence comparison, not to reproduce a high-coupling material. For quantitative predictions, refine the rotor harmonic rank and crystallite grid independently.
