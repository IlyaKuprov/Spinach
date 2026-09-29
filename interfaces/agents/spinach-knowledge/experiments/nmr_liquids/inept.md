# experiments/nmr_liquids/inept.m

- Signature: `fid=inept(spin_system,parameters,H,R,K)`

## Purpose

This is the non-refocused INEPT variant: the source says it returns the directly acquired coupled antiphase spectrum, rather than a refocused, broadband-decoupled variant. It cites [this paper](https://doi.org/10.1021/ja00497a058).

## Sequence and signal

The source starts from isotropic thermal equilibrium and detects the first configured nucleus. It applies a 90-degree x pulse to the second nucleus, evolves for `tau=abs(1/(4*parameters.J))`, applies simultaneous 180-degree y pulses to both configured nuclei, evolves for a second `tau`, then applies a 90-degree x pulse to the first nucleus. Opposite 90-degree y phases on the second nucleus are combined by a difference phase cycle before direct acquisition. The code does not specify an explicit coherence-order filter, so no narrower pathway selection is asserted here. With `J` in Hz, `tau` is in seconds; dwell time is `1/parameters.sweep` for a sweep width in Hz. The output `fid` is a one-dimensional FID with `parameters.npoints` samples. The source notes that `dilute.m` can generate carbon isotopomers.

## Inputs

- `parameters.sweep`: one positive sweep width in Hz.
- `parameters.npoints`: positive integer number of FID points.
- `parameters.spins`: two different working isotope names in a cell array, ordered `{F1 F2}`; source example: `{'15N','1H'}`.
- `parameters.J`: working scalar coupling in Hz.
- `H`, `R`, `K`: Hamiltonian matrix, relaxation superoperator, and kinetics superoperator supplied by the context function. The source requires the `sphten-liouv` formalism and same-sized matrices.

[Canonical MATLAB source](https://github.com/IlyaKuprov/Spinach/blob/main/experiments/nmr_liquids/inept.m) · [Spin Dynamics Wiki page](https://spindynamics.org/wiki/index.php?title=inept.m).
