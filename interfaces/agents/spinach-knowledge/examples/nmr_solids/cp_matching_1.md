# examples/nmr_solids/cp_matching_1.m

Signature: cp_matching_1()

## What the example models

The file describes a Hartmann–Hahn cross-polarisation (CP) matching test between ¹H and ¹⁵N under MAS. Its explicit system has two spins, ¹H and ¹⁵N, sys.magnet=9.394, scalar Zeeman entries 0.1495 and 0, and coordinates [−1.11551509, 1.65289357, −1.19927242] and [−2.67552180, 0.95825426, 0]. The basis is sphten-liouv with approximation none. Coordinates and Zeeman values are code settings; the file does not spell out additional interaction terms or their units.

## MAS and CP scan

The example passes rate 10000, rotor-axis vector [sqrt(2/3), 0, sqrt(1/3)], max_rank 3, and powder grid rep_2ang_200pts_oct to the experiment call. It starts from the ¹H Lx state, detects with the ¹⁵N Lx coil, sets zero excitation operators, and uses ten time steps each set to 4e-5. The contact experiment defines the time-step values in seconds, giving 400 µs total contact time; the example does not state units for rate or the second irradiation setting.

It samples 120 ¹H irradiation settings from 20e3 to 80e3, holding the ¹⁵N irradiation setting at 50e3. The scan runs in a parfor loop. The plot converts the scanned values by 1e3 and labels the horizontal axis as ¹H spin-lock RF power in kHz, so the displayed range is 20–80 kHz. At each setting the code calls singlerot with cp_contact_hard, stores real(fid(end)), and plots that final contact-curve ¹⁵N signal in a.u. against the scanned ¹H power.

## Scope of the result

This is a computed final-contact signal-versus-power curve, not an experimental spectrum or a validation result. `cp_contact_hard` returns the coil expectation before contact and after each of ten 40 µs slices: eleven points per setting, with `fid(end)` at 400 µs. No further sequence internals are inferred here. The source header estimates calculation time as seconds; that is a source comment, not a timed run in this note.

Source: https://github.com/IlyaKuprov/Spinach/blob/main/examples/nmr_solids/cp_matching_1.m
