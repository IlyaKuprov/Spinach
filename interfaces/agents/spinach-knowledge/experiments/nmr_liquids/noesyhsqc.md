# experiments/nmr_liquids/noesyhsqc.m

[MATLAB source](https://github.com/IlyaKuprov/Spinach/blob/main/experiments/nmr_liquids/noesyhsqc.m) · [Spin Dynamics Wiki](https://spindynamics.org/wiki/index.php?title=noesyhsqc.m)

- Signature: `fid=noesyhsqc(spin_system,parameters,H,R,K)`

## Purpose and sequence

A phase-sensitive NOESY-HSQC sequence hard-wired to `{F1,F2,F3}={1H,15N,1H}`, with carbon decoupled throughout when `13C` is present. The source cites [the NOESY-HSQC reference](https://doi.org/10.1021/bi00441a004), [the earlier precursor](https://doi.org/10.1016/0022-2364(90)90227-Z), [the handbook chapter](https://doi.org/10.1002/9780470034590.emrstm0563), and [the later article](https://doi.org/10.1016/j.jmr.2014.04.002). It initialises proton longitudinal magnetisation internally and uses proton L+ detection; the header notes that the relaxation destination is the zero state, not a thermalised state.

The code evolves the NOESY t1 period with nitrogen decoupled, explicitly branches the proton coherence into +1 and −1 components, applies 90° x pulses, homospoils, and evolves the NOESY mixing time under decoupled relaxation and kinetics. It then performs the HSQC transfer with two periods of `delta=1/(4*abs(J))`, combined proton/nitrogen 180° pulses, and the listed y/x pulses, before applying the final nitrogen 90° pulse. In the backward detection trajectory it selects zero proton coherence and splits nitrogen coherence into +1 and −1. Four forward/backward combinations are stitched into the States-processing fields; this description follows the explicit coherence projections and operations in the source rather than assigning additional pathways.

## Inputs

- `parameters.npoints`: three positive integer point counts ordered [t1 t2 t3].
- `parameters.sweep`: three positive sweep widths in Hz ordered [f1 f2 f3].
- `parameters.J`: HSQC-stage J coupling in Hz; the transfer delay is `abs(1/(4*J))` seconds.
- `parameters.tmix`: NOESY mixing time in seconds.
- `H`, `R`, and `K`: Hamiltonian, relaxation, and kinetics matrices from the context function. The required isotope set is `1H` and `15N`; `13C` is decoupled if present.

## Output

The four fields `fid.pos_pos`, `fid.pos_neg`, `fid.neg_pos`, and `fid.neg_neg` are used for subsequent States quadrature processing. Each sampled three-dimensional FID is permuted by the source to [t3 t2 t1], i.e. array dimensions `[npoints(3), npoints(2), npoints(1)]`.
