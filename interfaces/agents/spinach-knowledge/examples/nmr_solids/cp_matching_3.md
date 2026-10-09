# examples/nmr_solids/cp_matching_3.m

Signature: cp_matching_3()

## What the example computes

This example makes a two-parameter CP matching map for a ¹H–¹⁵N pair under MAS; it is a grid of RF settings, not a two-dimensional acquired NMR spectrum. The source sets sys.magnet=9.394, isotopes ¹H and ¹⁵N, Zeeman scalar entries 0.1495 and 0, and coordinates [−1.11551509, 1.65289357, −1.19927242] and [−2.67552180, 0.95825426, 0]. It uses the sphten-liouv basis with approximation none. No further interaction terms or units for these values are stated in this file.

## Rotor, powder and RF settings

The passed experiment parameters include rate 10000, axis [sqrt(2/3), 0, sqrt(1/3)], max_rank 3, powder grid rep_2ang_200pts_oct, ¹H Lx initial state, ¹⁵N Lx detection coil, zero excitation operators, and ten time-step entries of 4e-5. The contact experiment defines time-step durations in seconds, giving 400 µs total contact time; the example does not state rate units. Each RF axis uses 50 values from 0e3 to 50e3: the first row of irr_powers varies ¹H and the second varies ¹⁵N. Both plot labels use Hz, so each code range is 0–50,000 Hz.

For each ¹H value, a parallel inner sweep evaluates all ¹⁵N values through singlerot with cp_contact_hard. The matrix stores real(fid(end)) as cp(n,k), where n is the first (¹H) setting and k the second (¹⁵N) setting. The image is rendered with imagesc and the source labels its horizontal axis ¹H spin-lock RF power, Hz, and vertical axis ¹⁵N spin-lock RF power, Hz. This note preserves both the matrix indexing and labels as written; it does not infer an axis correction or claim a validated interpretation. The output is a simulated final-contact signal map, not a measured spectrum. `cp_contact_hard` returns the initial coil expectation followed by ten 40 µs contact samples; `fid(end)` is the eleventh point at 400 µs. The source header estimates calculation time as hours; this is not a timing measurement made here.

Source: https://github.com/IlyaKuprov/Spinach/blob/main/examples/nmr_solids/cp_matching_3.m
