# kernel/kinetics/kinetics.m

- Signature: `K=kinetics(spin_system)`
- Direct MATLAB source: [`kernel/kinetics/kinetics.m`](https://github.com/IlyaKuprov/Spinach/blob/main/kernel/kinetics/kinetics.m)
- Existing Wiki: [`kinetics.m`](https://spindynamics.org/wiki/index.php?title=kinetics.m)

## Purpose and output

Builds the chemical-kinetics superoperator from `spin_system.chem`. `K` is initialised with `mprealloc(spin_system,0)` and has the system superoperator dimensions. In a manually assembled Liouvillian, the source header specifies `1i*K`, for example `L=H+1i*R+1i*K`; Spinach context functions include kinetics automatically.

## Chemical-reaction block mapping

For each nonzero entry of `spin_system.chem.rates`, the source obtains source and destination parts from `spin_system.chem.parts`. It makes masks for basis rows involving those spins and requires the selected source and destination basis submatrices to be identical. With those row-index sets `S` and `D`, the exact insertion is `K[S,D] += rate * ones(|S|,|D|)`. Reaction rates are used directly; this function applies no unit conversion. They must use the inverse-time convention compatible with the propagator to which `K` is added.

## Magnetisation-flux mapping

Each nonzero `flux_rate` entry identifies source and destination spins. The routine is available for this contribution only in `sphten-liouv` formalism. It separates basis rows with one active spin from rows with more than one, removes rows active on both the source and destination (stationary states), and checks that the remaining subspaces match.

For the single-spin rows, with source indices `S`, destination indices `D`, and flux rate `f`, the code adds `+f` to `K[D,S]` and `-f` to the source diagonal `K[S,S]`. For multi-spin rows, `intramolecular` flux applies the same mapped transfer and source loss, preserving correlations; `intermolecular` flux applies source-diagonal loss only, damping those correlation orders. Only `intramolecular` and `intermolecular` are accepted flux types.

## Radical-pair recombination

When `chem.rp_theory` is nonempty and not `off`, the source constructs singlet/triplet projector superoperators from the configured electron indices. The supported branches are `exponential`, `haberkorn`, and `jones-hore`: respectively, subtracting the summed rate times the identity; subtracting the half-weighted left-plus-right singlet/triplet projector terms; or subtracting the summed-rate identity and adding the rate-weighted triplet and singlet left/right products. This branch is guarded to `sphten-liouv` or `zeeman-liouv` formalism; an unknown model errors.

## Input and integrity guards

Chemical reaction and flux terms each explicitly require `sphten-liouv`. Reaction source/destination basis subspaces must match; flux subspaces must match after stationary rows are excluded. An unrecognised flux type errors. Radical-pair recombination has its separate formalism and model guards. The source uses configured rates as supplied and does not normalise them or convert between Hz and angular-frequency units; no such unit conversion is present in this routine.

## Implementation lifecycle

The routine initialises `K`, adds reaction blocks, processes fluxes if present, then adds the configured radical-pair term if enabled. If the assembled matrix has no nonzero entries, it reports that no significant kinetics was specified.
