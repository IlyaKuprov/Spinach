# examples/nmr_overtone/cpmas_glycine_accum.m

- Signature: `cpmas_glycine_accum()`
- Source: [examples/nmr_overtone/cpmas_glycine_accum.m](https://github.com/IlyaKuprov/Spinach/blob/main/examples/nmr_overtone/cpmas_glycine_accum.m)

## What this example models

This example calculates a simulated 14N-overtone/proton cross-polarisation accumulation profile for glycine under MAS as contact duration changes. Its source comment calls the powder grid rough and attributes the glycine quadrupolar-tensor data to O'Dell and Ratcliffe, DOI [10.1016/j.cplett.2011.08.030](https://doi.org/10.1016/j.cplett.2011.08.030). The wrapper does not load measured spectra.

The modeled system has isotopes 14N and 1H and a field of 14.10220742 T. The nitrogen quadrupolar input is passed as `eeqq2nqi(1.18e6,0.53,1,[0 0 0])`; the existing KB describes the first value as 1.18 MHz and the asymmetry as 0.53, while the wrapper itself gives the bare helper arguments. Scalar entries are `{32.4 0}` and coordinates are `[0 0 0]` and `[0 0 1]`. The previous KB described the latter as 1 A apart; this wrapper lists coordinate values but does not label their units. It likewise does not state units for the scalar entries.

Relaxation is damped, with diagonal retention, zero equilibrium, and `damp_rate=300` (no unit is annotated here). The basis is `sphten-liouv` with no approximation. The code disables Krylov and trajectory-level options. Its grid identifier is `rep_2ang_200pts_oct`, and the source comment describes the grid as a rough powder grid.

## MAS, RF preparation/detection, and acquisition axis

The wrapper sets `theta=atan(sqrt(2))` and the rotor axis to `[sqrt(2/3) 0 sqrt(1/3)]`. It sets `rate=-19840`, max rank 7, sweep `[4.4e4 5.2e4]`, 256 acquired and zero-filled points, averaging method, and `axis_units='kHz'`. The plotting window is 44 to 52 kHz. The wrapper does not annotate units for `rate` or the sweep vector itself.

The proton initial state and 14N coil state are each built from the corresponding `Lz` and `Lx` terms weighted by the magic-angle sine and cosine; the wrapper also defines matching 14N and proton operators. RF settings are `rf_pwr=2*pi*[55.0e3 35.1e3]/sin(theta)` and `rf_frq=48e3`; their units are not documented at this call site. For each of ten iterations, `rf_dur=1e-5*n`, giving contact durations from 10 to 100 microseconds. Each iteration calls `singlerot(spin_system,@overtone_cp,parameters,'qnmr')` and plots the real spectrum in one of ten panels.

## Scope and limits

The wrapper exposes the prepared/detected operators, MAS settings, RF parameters, contact-duration sweep, simulation call, and plotted spectra. It does not show the internal pulse program behind `@overtone_cp`, so pulse ordering, gradients, and receiver cycling are not inferred. The source comment estimates minutes of computation; no measured runtime or experimental comparison is reported.
