# examples/nmr_liquids/noe_four_spin.m

- Signature: `noe_four_spin()`
- Source: [`examples/nmr_liquids/noe_four_spin.m`](https://github.com/IlyaKuprov/Spinach/blob/main/examples/nmr_liquids/noe_four_spin.m)

## Purpose and model

An inversion-recovery NOE-effect simulation on a simple four-proton spin system. The source describes inversion of the rightmost proton, followed by pulse-acquire detection after a five-second mixing period; it characterises sequential NOE hops with alternating signs as the intended visible effect. This is a source description, not an independently supplied measured spectrum. The calculation-time estimate is seconds.

The four chemical-shift entries are `[1.0, 2.0, 3.0, 4.0]`; nearest-neighbour scalar couplings are set to 10 for pairs (1,2), (2,3), and (3,4). Coordinates place the spins at z = 0, 2, 4, and 6 with x = y = 0. The field parameter is 14.1. The source does not annotate units for the shifts, couplings, coordinates, or field parameter. The relaxation model is Redfield with Di Bari equilibrium, `rlx_keep='kite'`, temperature parameter 298, and correlation-time entry `200e-12` (units are not annotated).

## Inversion, mixing, and acquisition

The code builds the relaxation superoperator and thermal-equilibrium state, then constructs an inverted-spin initial state using `Lz` for spin 1. It evolves for 5.0 s under the relaxation superoperator and subtracts the unperturbed equilibrium state to isolate the NOE deviation. The subsequent pulse-acquire experiment detects proton `L+`, uses a `Ly` pulse of `pi/2`, and applies no decoupling. Acquisition parameters are offset 1400, sweep 4500, 8192 points, zero-filling to 65536, ppm axis units, and axis inversion; the source does not state offset or sweep units.

The liquid simulation uses `@hp_acquire`; exponential apodisation with parameter 6 precedes a shifted FFT, and the real spectrum is plotted over 0.8–4.2 ppm. The source provides no experimental NOE rates, calibrated intensities, or DOI.
