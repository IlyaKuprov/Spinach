# cce_sic_divacancy

`[answer,diagnostics]=cce_sic_divacancy()` reconstructs the CCE-1 Hahn-echo traces in Supplementary Figure 5 of Seo et al., [doi:10.1038/ncomms12935](https://doi.org/10.1038/ncomms12935), for the axial (kk) divacancy in natural 4H–SiC. It returns normalised coherence with time along the rows and magnetic field along the columns. Diagnostics gives the time axis in seconds, field axis in tesla, and the lattice and distinct-site counts.

The `E3` electron represents the triplet, and the echo swaps the projected `m_s=+1` and `m_s=0` manifolds. The axial electron splitting contributes only scalar phases in these manifolds, which cancel in the ideal echo. The model retains longitudinal and transverse nuclear hyperfine components from a point electron dipole at the silicon vacancy and uses a 50 Å bath, with occupancies 4.7% for `29Si` and 1.1% for `13C`. The bulk hexagonal coordinates are from the axial-divacancy PyCCE tutorial; there are no relaxed defect coordinates or fitted contact terms.

For each lattice site, the empty-site factor is one and the occupied-site factor is the projected Spinach `clexpn` Hahn echo. At CCE-1 the factors multiply, so independent Bernoulli substitution gives the isotope average directly as a product of single-site averages, without drawing a finite set of random baths. At fixed isotope, basal distance and absolute axial height determine the point-dipole singleton response in a c-axis field; equivalent sites share one calculation.

The echo modulations are electron spin echo envelope modulation (ESEEM), not a decay envelope from nuclear pair flip-flops. Neither the CCE-2 traces in the main article nor their fitted coherence times are targets of this example.
