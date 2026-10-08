# examples/nmr_liquids/hmbc_camphor.m

- MATLAB implementation: [examples/nmr_liquids/hmbc_camphor.m](https://github.com/IlyaKuprov/Spinach/blob/main/examples/nmr_liquids/hmbc_camphor.m)

- Signature: `hmbc_camphor()`

## Purpose and spin system

Liquid-state HMBC for camphor at natural 13C abundance; the source comment estimates seconds of calculation time. The wrapper parses `../standard_systems/camphor.log` and maps H/C to `1H`/`13C` through `g2spinach`; its comment identifies the coordinates, shielding anisotropies, and J-couplings as vacuum-DFT-derived. It sets `options.min_j=3.0` and `options.no_xyz=0`, and passes the source range expression `[31.8-0.35 182.1+7.14]` to the importer.

## Acquisition and simulation

The wrapper sets `sys.magnet=14.1`, enables `zte` and `greedy` with proximity cutoff 4.0, and uses `sphten-liouv` / `IK-2`, scalar-coupling connectivity, and proximity level 1. It sets `J=140` Hz, `delta_b=60e-3` s, sweeps `[40000 1500]` Hz, offsets `[18000 900]` (units not stated here), `[128 128]` points, and `[512 512]` zero-fill points. `spins={'13C','1H'}` assigns carbon to F1 and proton to F2; the example uses Hz axes. `dilute(...,'13C')` supplies natural-abundance carbon isotopomers.

The wrapper calls `liquid(...,@hmbc,parameters,'nmr')`; the separate `experiments/nmr_liquids/hmbc.m` pulse program uses proton excitation and detection, carbon pulses, a J-set delay and `delta_b`, and indirect carbon evolution with proton decoupling. The wrapper supplies settings, not the internal pulse operations. Hamiltonian, relaxation, and kinetics superoperators are passed from the `liquid` context; this wrapper specifies no relaxation rates or model. Each isotopomer is simulated in `parfor`, cosine-apodised in both dimensions, Fourier-transformed and summed; the absolute spectrum is plotted in positive mode. The pulse-program references are https://doi.org/10.1021/ja00268a061 and https://doi.org/10.1016/0022-2364(88)90172-2.
