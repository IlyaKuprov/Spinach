# experiments/nmr_liquids/noesyhsqc.m

- Signature: `fid=noesyhsqc(spin_system,parameters,H,R,K)`

## Purpose

Phase-sensitive NOESY-HSQC. The sequence is hard-wired to {F1,F2,F3} = {1H,15N,1H}, with carbon decoupled throughout. The source cites [10.1021/bi00441a004](https://doi.org/10.1021/bi00441a004), [10.1016/0022-2364(90)90227-Z](https://doi.org/10.1016/0022-2364(90)90227-Z), [10.1002/9780470034590.emrstm0563](https://doi.org/10.1002/9780470034590.emrstm0563), and [10.1016/j.jmr.2014.04.002](https://doi.org/10.1016/j.jmr.2014.04.002). The first citation is the NOESY-HMQC precursor noted by the source.

## Sequence and output

The routine starts from proton longitudinal magnetization (`Lz`) and assumes an unthermalised relaxation superoperator whose destination is the zero state. It decouples carbon when present, discretizes three evolution dimensions, and uses a J-transfer interval `delta = abs(1/(4*parameters.J))`. The NOESY stage separates the two proton coherence pathways, applies homospoil and mixing under relaxation and kinetics, then performs HSQC transfer. Forward density trajectories and backward detection trajectories are joined with `stitch`; the resulting dimensions are permuted to [F3 F2 F1].

- Output: `fid`, with fields `pos_pos`, `pos_neg`, `neg_pos`, and `neg_neg` for subsequent States quadrature processing.
- The implementation requires `sphten-liouv` formalism and both 1H and 15N isotopes in the spin system.

## Inputs

- `parameters.npoints`: three positive integer point counts ordered [t1 t2 t3].
- `parameters.sweep`: three positive sweep widths ordered [f1 f2 f3], in Hz.
- `parameters.J`: non-zero real HSQC-stage J-coupling in Hz.
- `parameters.tmix`: non-negative NOESY mixing time in seconds.
- `H`, `R`, and `K`: Hamiltonian, relaxation, and kinetics matrices supplied by the context function; they must have matching dimensions.

## Reference link

[Spinach Wiki: noesyhsqc.m](https://spindynamics.org/wiki/index.php?title=noesyhsqc.m)
