# experiments/nmr_liquids/ct_hsqc.m

[MATLAB source](https://github.com/IlyaKuprov/Spinach/blob/main/experiments/nmr_liquids/ct_hsqc.m) · [Spinach Wiki: ct_hsqc.m](https://spindynamics.org/wiki/index.php?title=ct_hsqc.m)

- Signature: `fid=ct_hsqc(spin_system,parameters,H,R,K)`

## Purpose

Constant-time phase-sensitive HSQC, citing [DOI 10.1016/0022-2364(92)90144-V](https://doi.org/10.1016/0022-2364(92)90144-V) and [DOI 10.1007/BF00227470](https://doi.org/10.1007/BF00227470). The source describes a two-spin liquid-NMR sequence, not an imaging or flow experiment.

## Inputs

- `parameters.sweep`: two positive sweep widths in Hz, ordered [F1, F2].
- `parameters.npoints`: two point counts, ordered [F1, F2].
- `parameters.spins`: two different isotope names in a cell array, ordered {F1, F2}; the header's example is `{'13C','1H'}`.
- `parameters.J`: nonzero working scalar coupling in Hz. The source uses `abs(1/(4*parameters.J))` for its J-evolution interval.
- `parameters.decouple_f2`: optional cell array of isotope names to decouple in F2; if omitted, the source sets it to an empty cell array. The source header gives `{'15N','13C'}` as an example.
- `parameters.rho0`: optional initial state; if omitted, the source uses `Lz` on spin 2.
- `parameters.coil`: optional detection state; if omitted, the source uses `L+` on spin 2.
- `H`, `R`, and `K`: Hamiltonian, relaxation, and kinetics matrices supplied by the context function; the source combines them as `H + 1i*R + 1i*K`. The implementation requires the `sphten-liouv` formalism.

## Sequence, phase selection, and detection

The source's pulse operators act on spin 1 (named `Cx`) and spin 2 (`Hx`/`Hy`). Preparation starts from spin 2, applies a 90-degree x pulse, a J interval, simultaneous pi inversion on both spins, a second J interval, and a 90-degree y pulse on spin 2. A difference of positive and negative 90-degree x rotations on spin 1 establishes the phase-sensitive pathway. The F1 constant-time loop uses a time grid from `parameters.sweep(1)`, with delays `(CT-t1)/2`, `CT/2`, and `t1/2`, separated by the source's spin-1 and spin-2 pi pulses.

For each F1 state, the source separates the + and - pathways by selecting zero coherence on spin 2 with +1 or -1 coherence on spin 1. It applies the subsequent pulses and J evolution, decouples the requested F2 isotopes, and detects on spin 2 while evolving F2 at timestep `1/parameters.sweep(2)`. The source uses isotope dilution for natural-abundance simulations; see `dilute.m`.

## Output and scope

- `fid.pos` and `fid.neg`: the two components of the States quadrature signal, with the F1 states followed by F2 observable evolution.

This documents the parameterised pulse sequence; it does not report measured data or assert that a simulation was executed or validated.
