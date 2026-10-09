# examples/nmr_solids/cp_matching_2.m

Signature: cp_matching_2()

## What the example models

This is a Hartmann–Hahn ¹H–¹⁵N cross-polarisation matching test without exchange, described by the source as running with low power on ¹⁵N and showing matching-condition reflections with opposite phase. The code sets a two-spin system at sys.magnet=9.394, Zeeman scalar entries 0.1495 and 0, and coordinates [−1.11551509, 1.65289357, −1.19927242] and [−2.67552180, 0.95825426, 0]. It uses the sphten-liouv basis with approximation none. The source supplies coordinates but does not define additional interaction terms or state their units.

## MAS and power sweep

The model passes rate 10000, axis [sqrt(2/3), 0, sqrt(1/3)], max_rank 3, and grid rep_2ang_200pts_oct. The initial state is ¹H Lx, the detection coil is ¹⁵N Lx, excitation operators are zero, and the time grid contains ten entries of 4e-5. The contact experiment defines the time-step values in seconds, giving 400 µs total contact time; the example does not state units for rate or the fixed ¹⁵N irradiation value.

A parfor loop scans fifty ¹H irradiation settings from 0e3 to 30e3; ¹⁵N is held at 1e3. The plotted x-axis divides the scanned values by 1e3 and is labelled ¹H spin-lock RF power in kHz, giving a displayed range of 0–30 kHz. For each setting, singlerot is called with cp_contact_hard; the stored observable is real(fid(end)), plotted as ¹⁵N signal in a.u. The opposite-phase reflection description is the source comment about the intended/illustrated pattern, not an independently measured or validated spectrum.

## Scope and call boundary

`cp_contact_hard` returns the coil expectation before contact and after each of ten 40 µs slices: eleven points per setting, with `fid(end)` at 400 µs. This is a simulated final-contact signal-versus-power sweep, not a full experimental sequence or a 2D acquired spectrum. The source header estimates calculation time as seconds; no timing measurement is reported here.

Source: https://github.com/IlyaKuprov/Spinach/blob/main/examples/nmr_solids/cp_matching_2.m
