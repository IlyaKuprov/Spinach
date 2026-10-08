# examples/parahydrogen/case_studies/hyperpolarised_deuterium/just_bubbling.m

- Signature: just_bubbling()
- Source: [examples/parahydrogen/case_studies/hyperpolarised_deuterium/just_bubbling.m](https://github.com/IlyaKuprov/Spinach/blob/main/examples/parahydrogen/case_studies/hyperpolarised_deuterium/just_bubbling.m)

## What is modelled

This companion calculation follows spin-state coefficients during ortho-D2 bubbling in a parahydrogenation-catalyst model, without an RF pulse or acquisition. Each deuteron has spin 1, so the two-spin basis includes singlet (S), triplet (T), and quintet (Q) sectors. Exchange between the D2 pair (spins 1–2) and the catalyst-associated pair (spins 3–4) is combined with quadrupolar evolution and relaxation; the model therefore asks how the imposed bubbling/pumping process redistributes these spin-state components, and how they evolve after bubbling stops. It does not report experimental state populations or explicitly simulate para-H2 spin transfer: the coded system contains deuterons only, and the effective pumping rate is a guess.

The field is 7.05 T; the shifts are [4.55, 4.55, −13.5, −16.5], and the scalar couplings are 12.0 Hz and 0.24 Hz within the respective pairs. The DFT NQI tensors for spins 3 and 4 are set to 10^3 times [[108.2, 0.1, 28.1], [0.1, −55.6, 4.5], [28.1, 4.5, −52.6]] and [[−55.2, 8.1, −14.5], [8.1, −6.3, −73.9], [−14.5, −73.9, 61.5]]. The code supplies positions for these two spins and none for the D2 pair. It defines two chemical compartments with concentrations [1, 0] and rate matrix [[−1, 5000], [1, −5000]]. The source does not state units beside the matrix or NQI entries; the factor 10^3 is retained here as coded rather than silently converting it.

## Evolution and plotted observables

The `sphten-liouv` model uses no basis approximation. Zero-equilibrium secular Redfield relaxation uses correlation times 1 ps and 400 ps, R1 inputs [0.04, 0.04, 0, 0], and R2 inputs [8, 8, 0, 0]. The initial density is the unit state. During bubbling, `magpump` targets the singlet and quintet components with its identity part removed and rate 0.1; the MATLAB comment says this rate is a guess needing a proper value. The resulting state then evolves under the free Hamiltonian, relaxation, and exchange. The plotted real projections are the singlet, selected triplet, and quintet coefficients, labelled in arbitrary units. This is a simulated trajectory, not an experimental hyperpolarisation measurement.

The trajectory calls use step intervals of 0.007 and 0.03 with 1,000 steps each, giving 7 seconds of bubbling followed by 30 seconds of free evolution, as also shown on the plotted 0–7 and 7–37 second axis. The plotted timing therefore reflects the integration step and step count; 0.007 and 0.03 are not the total segment durations. The MATLAB header says the paper link will follow; no DOI is provided.

Zero track elimination is explicitly enabled with `zte` in `sys.enable`.
