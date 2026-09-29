# examples/nmr_solids/cp_matching_4.m

Signature: cp_matching_4()

## What the example models

The source describes CP between ¹H and ¹⁵N with conformational exchange between two geometries whose N–H vectors differ by 90°. The code represents four spins as two separate ¹H–¹⁵N pairs, with the N positions at [0,0,2] and [0,2,0]; the two exchange groups are [1 2] and [3 4]. It sets sys.magnet=9.394 and four zero scalar Zeeman entries, uses exchange-rate matrix [−5000, +5000; +5000, −5000] and concentrations [1,1]. No units for these exchange-rate values are stated. The basis is sphten-liouv with approximation none.

## MAS and CP scan

The example passes rate 10000, axis [sqrt(2/3), 0, sqrt(1/3)], max_rank 3, and powder grid rep_2ang_200pts_oct. It uses a ¹H Lx initial state, ¹⁵N Lx detection coil, zero excitation operators, and ten time-step entries of 4e-5. The time-step and rate units are not stated in the source.

A parfor loop scans 60 ¹H irradiation settings from 20e3 to 80e3 while holding the ¹⁵N irradiation setting at 50e3. The plotted scan values are divided by 1e3 and the axis is labelled ¹H spin-lock RF power in kHz, giving a displayed range of 20–80 kHz. At each point, singlerot is called with cp_contact_hard; real(fid(end)) is plotted as ¹⁵N signal in a.u. The plotted quantity is a computed final FID sample, not a measured spectrum or a validation result. Contact duration and pulse/sequence details are delegated to the named call and are not specified here. The source header estimates calculation time as seconds, not a measured timing result.

Source: https://github.com/IlyaKuprov/Spinach/blob/main/examples/nmr_solids/cp_matching_4.m
