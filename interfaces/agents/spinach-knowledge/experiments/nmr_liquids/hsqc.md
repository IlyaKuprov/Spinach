# experiments/nmr_liquids/hsqc.m

- Signature: `fid=hsqc(spin_system,parameters,H,R,K)`

## Purpose

Phase-sensitive heteronuclear single-quantum coherence (HSQC), using the configured F1 and F2 nuclei. The source cites the [HSQC paper](https://doi.org/10.1016/0009-2614(80)80041-8) and [review](https://doi.org/10.1002/cmr.a.10095).

## Sequence and signal

The default initial state is longitudinal magnetisation on F2; the default receiver detects F2. The pulse sequence evolves for `abs(1/(2*parameters.J))` around its F1/F2 inversion pulses, applies phase-sensitive F2 and F1 pulses, and builds the indirect F1 evolution. It applies the requested midpoint F1 refocusing pulses, then selects the F1 coherence orders +1 and -1 while retaining F2 order 0. The two selected components receive the subsequent F1/F2 pulses and coupling delays, F2 decoupling, and F2 detection. This description follows the operations and coherence filter in the source; it does not assign a more specific spin-operator pathway.

The dwell times are `1./parameters.sweep(1)` in F1 and `1./parameters.sweep(2)` in F2, with sweep widths in Hz. The source returns the States quadrature components as `fid.pos` and `fid.neg`. Each is a two-dimensional array with `npoints(2)` F2 acquisition samples in rows and `npoints(1)` F1 increments in columns; the source passes the F1 state stack directly to `evolution` without transposing the result.

## Inputs

- `parameters.sweep`: two positive sweep widths, `[F1 F2]`, in Hz.
- `parameters.npoints`: two positive integer point counts, `[F1 F2]`.
- `parameters.spins`: two different active isotope names in a cell array, ordered `{F1 F2}`; source example: `{'13C','1H'}`.
- `parameters.decouple_f2`: required cell array of isotope names to decouple during F2; source examples: `{'15N','13C'}`.
- `parameters.decouple_f1`: required cell array of isotope names receiving midpoint 180-degree refocusing pulses in F1; source examples: `{'1H','13C'}`.
- `parameters.J`: working scalar coupling in Hz; its absolute value sets the coupling delay above.
- Optional `parameters.rho0` and `parameters.coil`: initial state and detection state; if omitted, the source constructs the F2 longitudinal state and F2 raising-operator receiver.
- `H`, `R`, `K`: Hamiltonian matrix, relaxation superoperator, and kinetics superoperator supplied by the context function. The source requires the `sphten-liouv` formalism and same-sized matrices.

Natural-abundance simulations should use isotope dilution; see `dilute.m`.

[Canonical MATLAB source](https://github.com/IlyaKuprov/Spinach/blob/main/experiments/nmr_liquids/hsqc.m) · [Spin Dynamics Wiki page](https://spindynamics.org/wiki/index.php?title=hsqc.m).
