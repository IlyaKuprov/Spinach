# examples/parahydrogen/case_studies/hyperpolarised_deuterium/kinetic_isotope_effect.m

- Signature: kinetic_isotope_effect()
- Source: [examples/parahydrogen/case_studies/hyperpolarised_deuterium/kinetic_isotope_effect.m](https://github.com/IlyaKuprov/Spinach/blob/main/examples/parahydrogen/case_studies/hyperpolarised_deuterium/kinetic_isotope_effect.m)

## Model and scope

This calculation uses a four-deuteron (spin-1) model for ortho-D2 bubbling in the presence of a parahydrogenation catalyst. Exchange couples the D2 pair (spins 1–2) to a catalyst-associated pair (spins 3–4); spin evolution also includes the pair's singlet, triplet, and quintet structure, DFT-supplied NQI tensors, and secular Redfield relaxation. The title names a kinetic isotope effect, but the script contains only this deuterium model: it does not set up a protium control, vary an isotope-dependent rate, or calculate an isotope-effect ratio. Its plots are simulated state coefficients and coherences, not measured populations or a measured isotope effect. The coded system contains deuterons only; it does not explicitly model para-H2 proton spin transfer.

Inputs are a 7.05 T field, shifts [4.55, 4.55, −13.5, −16.5], and within-pair scalar couplings 12.0 Hz and 0.24 Hz. The spin-3 and spin-4 NQI tensors are respectively 10^3 times [[108.2, 0.1, 28.1], [0.1, −55.6, 4.5], [28.1, 4.5, −52.6]] and [[−55.2, 8.1, −14.5], [8.1, −6.3, −73.9], [−14.5, −73.9, 61.5]]. The associated DFT coordinates are [−1.98, 0.45, −0.55] and [−0.25, 1.33, −1.65], while the D2 pair has no coordinates. Chemical compartments use concentrations [1, 0] and rate matrix [[−1, 5000], [1, −5000]]. The source does not state units beside that matrix or the NQI tensor entries, so the numeric forms and 10^3 factor are reported as coded. The basis is `sphten-liouv` without approximation; relaxation combines Redfield and T1/T2 settings with zero equilibrium and secular retention, correlation times 1 ps and 400 ps, R1 inputs [0.04, 0.04, 0, 0], and R2 inputs [8, 8, 0, 0].

## State trajectory and pulse-created coherences

The initial density is the unit state. A bubbling Liouvillian adds a `magpump` term targeting singlet and quintet components, with rate 0.1; the source comments that this bubbling rate is a guess that needs a proper value. After storing the bubbling and free-evolution trajectories, the script projects the real singlet, selected triplet, and quintet coefficients. It also applies a 45-degree deuterium rotation to the stored trajectory and plots three transition-coherence projections: T1 to T0, Q1 to Q0, and Q2 to Q1. These are spin-state observables in arbitrary units, not a pulse-acquire NMR spectrum.

The trajectory calls use step intervals of 0.007 and 0.03 with 1,000 steps each, giving 7 seconds of bubbling followed by 30 seconds of free evolution, consistent with the plotted 0–7 and 7–37 second axis. These values are step intervals, not total segment durations. The MATLAB header says a paper link will follow, and supplies no DOI.

Zero track elimination is explicitly enabled with `zte` in `sys.enable`.
