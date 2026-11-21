# examples/cce/gcce_sic_clock.m

`gcce_sic_clock()` plots the free-induction coherence magnitude of the basal kh divacancy in natural 4H–SiC at the zero-field clock transition. It reproduces the mixed-bath gCCE1 and gCCE2 traces of Fig. 4(c) in Onizhuk and Galli, [doi:10.1002/adts.202100254](https://doi.org/10.1002/adts.202100254), not the Monte Carlo–averaged trace in Fig. 4(a).

The observed coherence is between the spin-one `mS=0` state and the even superposition of `mS=±1` selected by positive transverse zero-field splitting. The nuclear density is a maximally mixed product, with no exterior mean field. The 40 Å seed-1 configuration contains 161 carbon-13 and 593 silicon-29 nuclei, with 1,963 pairs inside 8 Å. The labelled `sic_clock_kh_sites.csv` lists the isotope, position in Å, and hyperfine tensor for each site (inner-shell GIPAW, outer-shell point dipoles). `sic_clock_kh_pairs.csv` lists each pair in full-system numbering (the electron is spin 1) and its nuclear dipolar tensor. Both tables give coupling tensors in Hz. The hyperfine columns retain PyCCE's nuclear-first Cartesian index order, transposed on assignment to Spinach. No full-bath basis is formed.

The plot displays only the normalised gCCE1 and gCCE2 magnitudes, without clipping. This is one finite isotope configuration and does not establish a bath-averaged coherence time.

[MATLAB source](../../../../../examples/cce/gcce_sic_clock.m)
