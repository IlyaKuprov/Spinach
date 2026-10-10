# kernel/kinetics/kinetics.m

## Direct-sum chemistry

`K=kinetics(spin_system)` compiles `chem.reactions` into sparse drains and product-row maps. Numeric first-order records return a constant sparse matrix. Numeric zero-rate higher-order records also permit the static route; an entirely zero-rate network yields a sparse zero matrix usable by ordinary linear contexts. Rate callbacks remain dynamic even if a sampled value is zero. Higher-order reactions or time-dependent rate handles return `K(t,eta)`, evaluated on the instantaneous concentration-weighted state. Chemistry-free systems retain their zero generator in every formalism. Both Liouville formalisms support first-order and mass-action reactions; Hilbert reactions expose a first-order matrix derivative.

The dissipative convention is `L=H+1i*R+1i*K`. A state-dependent generator uses the existing `step` handle route, for example `{ @(t,eta)1i*K(t,eta), t, 'RKMK4' }`. No block is normalised or divided by concentration. Every reactant occurrence contributes a drain multiplied by the concentrations of the other occurrences. Repeated products contribute repeated fills. Spin-free reactants are dynamic pools, not fixed-concentration reservoirs.

## Arrival closures and selectors

The additive closure carries each reactant's internal spin orders and distributes the unit arrival equally over the reactant occurrences, so the reaction event is counted once. Cross-reactant orders are omitted. The product closure gathers source coordinates from all but the last occurrence into the sparse matrix acting on the last; it retains their polarisation products without constructing a tensor-product basis.

Named singlet/triplet selectors use the left/right electronic projectors: Haberkorn loss is half their sum, and arrival is their product followed by the matching map. Jones–Hore variants use identity minus the complementary projector product for the drain. User selector pairs are substance-local left/right product superoperators with the reactant Liouville block dimensions, including in `zeeman-hilb`.

`kinetics(spin_system,'report')` prints the network, closure, matched and traced spins. Compilation reports product rows missing a source descriptor. Maps and selector products are compiled once per call, not inside time stepping. For a space-times-spin state, the handle assembles an independent sparse chemistry block per voxel; transport is added separately. Each time-dependent rate is evaluated and validated once at the shared stage time, then reused across voxels.

## Zeeman representation and matrix actions

For Zeeman input, the same reaction compiler runs on a complete spherical-tensor basis for each substance. Its local maps are transformed with the physical `sphten2zeeman/D_n` matrix and its inverse. Numeric first-order maps are transformed once. A state-dependent map pulls back the instantaneous vector, evaluates the shared polynomial generator, and transforms its result independently per voxel. This small-system reference implementation pays for complete local bases and dense basis transformations; it does not duplicate the closure, matching, or selector physics.

With nonempty `zeeman-hilb` reactions, `K(t,rho)` returns the matrix derivative through `hilb_action`. It is not a Hamiltonian and must not be supplied as a conjugation generator to `step`. Numeric and time-dependent first-order rates are supported. More than one reactant occurrence raises `Spinach:kinetics:hilbertMassAction` with the complete message “mass-action chemistry is not supported in zeeman-hilb formalism.” A ket formalism remains storage-only and rejects reaction records.

## Retired mechanisms

Legacy rate, flux, and radical-pair fields are rejected by `create`; this routine only reads explicit reaction records. First-order exchange is a directed record with matching, untracked exponential loss has empty products, and radical-pair channels use selectors. A permutation reaction changes all mapped orders; it must not be assumed identical to a legacy phenomenological flux model that left some correlations stationary.

Nonempty wavefunction reaction records raise `Spinach:kinetics:wavefunction`, naming zeeman-wavef explicitly; this covers first-order and mass-action chemistry.

Legacy global `bas.basis` matrices and `bas.irrep` fields are rejected at this entry point with named errors pointing to per-substance `bas.basis{n}`/`bas.offsets` and `bas.sym_fact(n)` symmetry data. Compiled structures remain ordinary MATLAB structs; arbitrary external dot reads are not intercepted.
The kinetics consumer also rejects all six retired chemistry fields, even when empty, with `Spinach:kinetics:retiredChemistry` pointing to `chem.reactions`; this covers manually modified or saved legacy structures that bypass `create`.
Reporting accepts mixed row/column substance memberships by formatting local spin lists as rows; the stored memberships and generator are unchanged.
