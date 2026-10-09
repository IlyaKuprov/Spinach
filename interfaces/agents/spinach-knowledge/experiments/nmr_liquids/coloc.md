# experiments/nmr_liquids/coloc.m

- MATLAB source: [experiments/nmr_liquids/coloc.m](https://github.com/IlyaKuprov/Spinach/blob/main/experiments/nmr_liquids/coloc.m)
- Spinach Wiki: [coloc.m](https://spindynamics.org/wiki/index.php?title=coloc.m)
- Sequence reference: [DOI 10.1016/0022-2364(84)90136-7](https://doi.org/10.1016/0022-2364(84)90136-7)

## Purpose

COLOC as shown in Fig. 1b of the cited paper, with the dashed pulses during the delta2 interval omitted. The code prepares an embedded echo in the F1 evolution period, selects proton coherence, transfers it during delta2, then detects on the second spin. The output is an FID intended for magnitude-mode processing; this description is of the implemented sequence, not a measured outcome. All sequence delays and readout evolve under `L=H+1i*R+1i*K`.

## Inputs and parameters

Signature: fid=coloc(spin_system,parameters,H,R,K)

- parameters.sweep: positive [F1 F2] sweep widths in Hz.
- parameters.npoints: positive integer point counts [F1 F2].
- parameters.spins: two isotope labels {F1 F2}; the source gives {'13C','1H'} as an example. The sequence starts on spin 1 and detects on spin 2.
- parameters.delta2: required COLOC delay in seconds; the source gives 40e-3 seconds as a typical value.
- parameters.delta1: optional echo delay in seconds. If omitted, the code sets it to (npoints(1)-1)/(2*sweep(1)), half the maximum F1 evolution time. A supplied value must be at least that large.
- H, R, and K: same-size Hamiltonian, relaxation, and kinetics matrices from the context function; the routine requires sphten-liouv formalism.

## Sequence and readout

The initial state is Lz on the first (F1) spin; the coil state is L+ on the second (F2) spin. After a 90-degree x pulse on F1, each F1 increment is embedded in an echo: evolution for half the increment, simultaneous 180-degree x pulses on both spins, and the remaining part of delta1. The source labels its selection block “Proton coherence selection”; the actual projection requests coherence order -1 on `parameters.spins{1}` (the first configured spin). A 90-degree F1 x pulse and F2 y pulse precede delta2; the code then selects order +1 on F2, decouples F1, and records the F2 observable evolution at 1/sweep(2) spacing for npoints(2) points.

For natural-abundance simulations, the source recommends isotope dilution via dilute.m.

No MATLAB execution or experimental signal is claimed here.
