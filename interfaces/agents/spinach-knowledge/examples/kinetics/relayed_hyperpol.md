# examples/kinetics/relayed_hyperpol.m

Source: https://github.com/IlyaKuprov/Spinach/blob/main/examples/kinetics/relayed_hyperpol.m

- Signature: `relayed_hyperpol()`
- Source citation: Christopher Pötzl, Figure S7, https://doi.org/10.1016/j.jmr.2024.107727

## Model and exchange

This model follows relayed NOE from hyperpolarised water to an ALA-GLY dipeptide. It contains 30 1H spins: ten peptide protons with coordinates (four marked labile and six aliphatic) and 20 water protons with empty coordinates. The missing water coordinates intentionally remove direct intermolecular cross-relaxation in this model. The ten peptide shifts are [8.45, 8.45, 8.45, 8.11, 3.73, 0.99, 0.99, 0.99, 3.99, 3.32] ppm; water shifts are 4.5 ppm. The four labile peptide spins (1-4) and the water spins (11-20) are linked by symmetric intermolecular flux-rate blocks set to `20`; the source does not annotate units for this value.

At 16.4 T and 298 K, the source selects Redfield T1/T2 relaxation, Dibari equilibrium, secular retention, and a correlation time of 1.2e-10 s. Peptide R1/R2 entries are zero and water entries are 0.1 Hz. The basis uses sphten-liouv, IK-1, full-tensor connectivity, proximity level 3, and interaction level 1. Hamiltonian, relaxation, and kinetics terms are assembled as `L=H+1i*R+1i*K`.

## Preparation, evolution, and output

Starting from isotropic equilibrium, the script replaces the `Lz` component only for the ten exchange-coupled water protons (spins 11–20) with a fully polarised `Wz` component. The other ten water protons (21–30) retain equilibrium populations and have no exchange-flux entries. It detects the aliphatic methyl protons 6-8 and H-alpha proton 5, then calls multichannel `evolution` with dt=0.125 s and 128 steps (16 s total). The plot shows the real CH3 and H-alpha magnetisation traces in arbitrary units.

The modeled relay is constrained by the coordinate-free water pool and specified exchange matrix; this page does not infer an experimental outcome beyond the source's stated Figure S7 context.
