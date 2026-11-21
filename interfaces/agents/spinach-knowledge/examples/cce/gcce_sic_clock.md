# examples/cce/gcce_sic_clock.m

`[answer,diagnostics]=gcce_sic_clock()` calculates the free-induction coherence of the basal kh divacancy in natural 4H–SiC at the zero-field clock transition. It reproduces the mixed-bath gCCE1 and gCCE2 traces of Fig. 4(c) in Onizhuk and Galli, [doi:10.1002/adts.202100254](https://doi.org/10.1002/adts.202100254), not the Monte Carlo–averaged trace in Fig. 4(a).

The observed coherence is between the spin-one `mS=0` state and the even superposition of `mS=±1` selected by positive transverse zero-field splitting. The nuclear density is a maximally mixed product, with no exterior mean field. The 40 Å seed-1 configuration contains 161 carbon-13 and 593 silicon-29 nuclei, with 1,963 pairs inside 8 Å. The physical-only fixture stores the isotopes, positions in Å, hyperfine tensors (inner-shell GIPAW, outer-shell point dipoles), pair indices in full-system numbering, and nuclear dipolar tensors. All tensor entries are in Hz. No full-bath basis is formed.

`answer` is the column of normalised complex gCCE2 coherences. `diagnostics.parameters.times` gives the time axis in seconds; `diagnostics.reduction.by_order` stores the central-only, gCCE1, and gCCE2 coherences in order, with matching validity masks and absolute numerical error bounds. The plot displays the last two magnitudes, without clipping. This is one finite isotope configuration and does not establish a bath-averaged coherence time.

[MATLAB source](../../../../../examples/cce/gcce_sic_clock.m)
