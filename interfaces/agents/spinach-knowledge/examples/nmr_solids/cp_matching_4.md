# examples/nmr_solids/cp_matching_4.m

Signature: cp_matching_4()

## What the example models

The source describes CP between ¹H and ¹⁵N with conformational exchange between two geometries whose N–H vectors differ by 90°. The code represents four spins as two separate ¹H–¹⁵N pairs, with the N positions at [0,0,2] and [0,2,0]; the two exchange groups are [1 2] and [3 4]. It sets sys.magnet=9.394 and four zero scalar Zeeman entries, uses two directed records with population generator [−5000, +5000; +5000, −5000] and concentrations [1,1]. No units for these exchange-rate values are stated. The basis is sphten-liouv with approximation none.

## MAS and CP scan

The example passes rate 10000, axis [sqrt(2/3), 0, sqrt(1/3)], max_rank 3, and powder grid rep_2ang_200pts_oct. It uses a ¹H Lx initial state, ¹⁵N Lx detection coil, zero excitation operators, and ten time-step entries of 4e-5. The contact experiment defines the time steps in seconds, giving 400 µs total contact time; the example does not state rate units.

A parfor loop scans 60 ¹H irradiation settings from 20e3 to 80e3 while holding the ¹⁵N irradiation setting at 50e3. The plotted scan values are divided by 1e3 and the axis is labelled ¹H spin-lock RF power in kHz, giving a displayed range of 20–80 kHz. At each point, singlerot is called with cp_contact_hard; real(fid(end)) is plotted as ¹⁵N signal in a.u. The plotted quantity is the last of eleven contact-curve samples: `cp_contact_hard` returns the initial coil expectation and one sample after each of ten 40 µs slices, so `fid(end)` is at 400 µs. This is not a measured spectrum or a validation result. The source header estimates calculation time as seconds, not a measured timing result.

Source: https://github.com/IlyaKuprov/Spinach/blob/main/examples/nmr_solids/cp_matching_4.m

The rate matrix is represented by two explicit first-order reaction records. Atom matching pairs equal-position spin indices in the two declared parts in both directions; the permuted geometry/tensors represent the conformational exchange. Detection uses unweighted `coil_state`.
