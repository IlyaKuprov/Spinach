# examples/nmr_solids/cp_matching_4.m

- Signature: `cp_matching_4()`

## Purpose

Tests the ¹H–¹⁵N Hartmann–Hahn power match in a model with conformational exchange between two geometries whose N–H vectors differ by 90°. The source estimates a calculation time of seconds.

## Physical / mathematical content

The model contains two ¹H–¹⁵N pairs, one for each geometry, with the two spin pairs grouped into exchanging chemical states. The exchange-rate matrix has off-diagonal rates of 5,000 s⁻¹ and equal state concentrations. Under a 10 kHz MAS rate, the experiment starts from ¹H transverse magnetisation and observes the ¹⁵N signal.

## Numerical / algorithmic content

The source uses the full `sphten-liouv` basis and `singlerot` with `@cp_contact_hard`. At MAS axis `[sqrt(2/3) 0 sqrt(1/3)]`, it scans 60 proton spin-lock powers from 20 to 80 kHz with the ¹⁵N power fixed at 50 kHz, using the `rep_2ang_200pts_oct` grid, `max_rank=3`, and ten 40 μs steps. Each scan point runs in a `parfor` loop; the real final FID point is plotted against ¹H power.

## Implementation structure

- Defines the two orientations as separate ¹H–¹⁵N pairs and specifies the chemical-exchange groups, rates, and concentrations.
- Builds the basis and transverse operators, then configures MAS, initial state, detection coil, and time grid.
- Performs the parallel RF-power sweep and plots the final ¹⁵N signal.
