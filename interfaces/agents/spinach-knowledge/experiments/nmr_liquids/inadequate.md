# experiments/nmr_liquids/inadequate.m

- Signature: `fid=inadequate(spin_system,parameters,H,R,K)`

## Purpose

INADEQUATE selects double-quantum coherence from coupled carbon pairs and returns a free induction decay. The source notes that at natural-abundance 13C this yields only 13C pair subspectra and recommends `dilute.m` for generating carbon-pair isotopomers. The implementation cites [this paper](https://doi.org/10.1021/ja00534a056).

## Sequence and signal

The sequence starts with longitudinal magnetisation on the configured nucleus and detects that same nucleus. It uses two delays `tau=abs(1/(4*parameters.J))` around a 180-degree y pulse, followed by a 90-degree x pulse and an explicit filter retaining coherence orders +2 and -2 for the configured nucleus. A final 90-degree x pulse precedes detection. The delay is in seconds when `J` is supplied in Hz. Acquisition uses dwell time `1/parameters.sweep` and returns `fid`, a one-dimensional FID with `parameters.npoints` samples.

## Inputs

- `parameters.sweep`: one positive sweep width in Hz.
- `parameters.npoints`: positive integer number of FID points.
- `parameters.spins`: active nucleus in a cell array; source example: `{'13C'}`.
- `parameters.decouple`: required cell array of nuclei to decouple; source example: `{'1H'}`.
- `parameters.J`: working scalar coupling in Hz.
- `H`, `R`, `K`: Hamiltonian matrix, relaxation superoperator, and kinetics superoperator supplied by the context function. The source requires the `sphten-liouv` formalism and same-sized matrices.

[Canonical MATLAB source](https://github.com/IlyaKuprov/Spinach/blob/main/experiments/nmr_liquids/inadequate.m) · [Spin Dynamics Wiki page](https://spindynamics.org/wiki/index.php?title=inadequate.m).
