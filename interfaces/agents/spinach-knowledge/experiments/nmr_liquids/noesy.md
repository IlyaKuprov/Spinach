# experiments/nmr_liquids/noesy.m

[MATLAB source](https://github.com/IlyaKuprov/Spinach/blob/main/experiments/nmr_liquids/noesy.m) · [Spin Dynamics Wiki](https://spindynamics.org/wiki/index.php?title=noesy.m)

- Signature: `fid=noesy(spin_system,parameters,H,R,K)`

## Purpose and sequence

A phase-sensitive homonuclear NOESY experiment. The source cites [the early NOESY paper](https://doi.org/10.1063/1.438208), [this report](https://doi.org/10.1016/0006-291X(80)90695-6), and [the later paper](https://doi.org/10.1016/0022-2364(82)90279-7). The code uses a four-step pulse phase cycle: after the first 90° x pulse, the second pulse alternates x/y with signs +π/2,+π/2,−π/2,−π/2; the third pulse is y with +π/2 in each step. It evolves t1, performs the mixing period, then acquires t2 with an L+ detection state. The returned quadrature channels subtract the third and fourth cycle members from the first and second to eliminate axial peaks. No explicit coherence-order projection is performed in this function.

By default, a homospoil step destroys all but longitudinal magnetisation before mixing, and the mixing evolution uses relaxation and kinetics, `iR+iK`. Setting `parameters.oldschool` to 1 disables homospoil and instead evolves the mixing period under the full `H+iR+iK` generator. The source describes this implementation as laid out for low memory use in extreme protein and nucleic-acid simulations rather than CPU speed; it notes that analytical decoupling is meaningful only in `sphten-liouv` formalism.

## Inputs

- `parameters.sweep`: two positive sweep widths in Hz; `parameters.npoints`: two positive integer point counts, both ordered by dimension.
- `parameters.spins`: nuclei label cell array; header examples include `{'1H'}` and `{'13C'}`.
- `parameters.tmix`: mixing time in seconds.
- Optional `parameters.decouple`: labels such as `{'13C','1H'}` or spin indices such as `[1 2]`.
- Supply `parameters.rho0` as the initial state, or omit it and set `parameters.needs={'rho_eq'}` to start from exact thermal equilibrium.
- Optional `parameters.oldschool`: set to 1 to disable the default homospoil gradient.
- `H`, `R`, and `K`: Hamiltonian, relaxation, and kinetics from the context function.

## Output

`fid.cos` and `fid.sin` are the two F1 hypercomplex components. Each is a two-dimensional FID with array shape `[npoints(2), npoints(1)]`: the direct t2 samples are rows and the t1 trajectory samples are columns.
