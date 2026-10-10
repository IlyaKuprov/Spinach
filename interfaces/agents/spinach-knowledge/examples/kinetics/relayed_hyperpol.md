# examples/kinetics/relayed_hyperpol.m

Source: https://github.com/IlyaKuprov/Spinach/blob/main/examples/kinetics/relayed_hyperpol.m

- Signature: `relayed_hyperpol()`
- Source citation: Christopher Pötzl, Figure S7, https://doi.org/10.1016/j.jmr.2024.107727

## Model and exchange

This model follows relayed NOE from hyperpolarised water to an ALA-GLY dipeptide. It contains 30 1H spins: ten peptide protons with coordinates (four marked labile and six aliphatic) and 20 water protons with empty coordinates. The missing water coordinates intentionally remove direct intermolecular cross-relaxation in this model. The ten peptide shifts are [8.45, 8.45, 8.45, 8.11, 3.73, 0.99, 0.99, 0.99, 3.99, 3.32] ppm; water shifts are 4.5 ppm. The peptide is one substance and each water proton is a separate unit-concentration pool. Forty additive A+B→A+B replacement records exchange each labile peptide spin (1–4) with each coupled water spin (11–20), at rate 20 and invariant unit concentrations. Matching retains the other peptide spins; tracing the departing spin destroys its correlations. Cross-molecule orders are omitted, as in the legacy intermolecular flux model.

At 16.4 T and 298 K, the source selects Redfield T1/T2 relaxation, Dibari equilibrium, secular retention, and a correlation time of 1.2e-10 s. Peptide R1/R2 entries are zero and water entries are 0.1 Hz. The peptide basis uses sphten-liouv, IK-1, full-tensor connectivity, proximity level 3, and interaction level 1; each water pool has a complete one-spin basis. All pools retain the same correlation time. The additive generator is frozen only because every concentration is invariant. Hamiltonian, relaxation, and kinetics terms are assembled as `L=H+1i*R+1i*K`.

## Preparation, evolution, and output

Starting from isotropic equilibrium, the script replaces the `Lz` component only for the ten exchange-coupled water protons (spins 11–20) with a fully polarised `Wz` component. The other ten water protons (21–30) retain equilibrium populations and have no replacement reactions. Unweighted `coil_state` vectors prepare the water projection and detect the aliphatic methyl protons 6-8 and H-alpha proton 5, then calls multichannel `evolution` with dt=0.125 s and 128 steps (16 s total). The plot shows the real CH3 and H-alpha magnetisation traces in arbitrary units.

The modelled relay is constrained by the coordinate-free water pool and specified replacement network; this page does not infer an experimental outcome beyond the source's stated Figure S7 context.

## Numerical comparison boundary

Mapped non-unit H and K agree exactly with the former flux representation. The stock global-unit DiBari correction introduces cross-substance non-unit relaxation entries even though unthermalised relaxation and left Hamiltonians have zero cross blocks. Per-substance relaxation deliberately omits that global-unit-mediated thermal coupling. Restoring only those stock cross-substance R blocks reduces the tight-propagation observable residual from `3.10306e-7` to `1.70802e-9`; this measures the dominant representation difference, not a defect in the published model.

The remaining residual is unresolved, so this example has no accepted `1e-10` parity verdict. A full-generator control rebuilt both propagators with identical scaling exponent 14 and `prop_chop=eps`; the cross-restored residual was `2.68137e-8`, not roundoff, while that stock control differed from its prior reduced propagation by `1.23923e-9`. Common scaling alone therefore did not establish parity. The open question is the remaining augmented unit-coordinate and propagation contribution after the measured cross-substance thermal term is isolated.
