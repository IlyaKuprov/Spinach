# experiments/nmr_liquids/cosy.m

- MATLAB source: [experiments/nmr_liquids/cosy.m](https://github.com/IlyaKuprov/Spinach/blob/main/experiments/nmr_liquids/cosy.m)
- Spinach Wiki: [cosy.m](https://spindynamics.org/wiki/index.php?title=cosy.m)
- Sequence references: [DOI 10.1063/1.432450](https://doi.org/10.1063/1.432450); [DOI 10.1016/0022-2364(82)90279-7](https://doi.org/10.1016/0022-2364(82)90279-7)

## Purpose

Phase-sensitive COSY with one analytically retained F1 coherence pathway. The implementation returns a two-dimensional free-induction decay. It describes a parameterised sequence, not an experimental or run-verified result. The sequence propagates states under `L=H+1i*R+1i*K`.

## Inputs and parameters

Signature: fid=cosy(spin_system,parameters,H,R,K)

- parameters.sweep: one positive sweep width in Hz, used for both dimensions; the sampling interval is 1/sweep seconds.
- parameters.npoints: two positive integer point counts [F1 F2].
- parameters.spins: one nucleus label in a cell array, e.g. {'1H'} or {'13C'}.
- parameters.angle: finite second-pulse angle in radians; pi/2 gives the usual 90-degree pulse, and the source also cites COSY45 and COSY60.
- H, R, and K: same-size Hamiltonian, relaxation, and kinetics matrices from the context function; the routine requires sphten-liouv formalism.

## Evolution, coherence selection, and detection

The routine starts from Lz magnetisation and an L+ detection state on the selected spin. It applies a 90-degree x pulse, records the F1 trajectory at 1/sweep spacing, and explicitly selects F1 coherence order +1. The second x pulse uses parameters.angle; direct F2 evolution is then detected with the L+ observable at the same reciprocal-sweep spacing. Thus the code retains a single phase-sensitive pathway rather than summing all F1 coherence orders.

The source advises magnitude-mode plotting when the second pulse differs from 90 degrees. No MATLAB execution or experimental signal is claimed here.
