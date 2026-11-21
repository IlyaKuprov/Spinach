# cce_hbn_quadrupoles

Syntax: `[answer,diagnostics]=cce_hbn_quadrupoles()`. The numeric answer is the normalised complex projected Hahn echo; diagnostics is the separate `clexpn` audit record.

The two-nucleus model uses explicit on-site quadrupole tensors for spin-3/2 `11B` and spin-1 `14N`, and a fixed magnetic-quantum-number transition. The complete order-two graph makes it possible to check the cluster expansion against independent dense bath propagators in `test_cce_examples`.

This model illustrates the quadrupole input convention, not the material-specific echo of the spin-1 boron vacancy. Haykal et al., [Fig. 4a,b](https://doi.org/10.1038/s41467-022-31743-0), report bulk gCCE3 Hahn-echo fits of 92 ns for the B-11-enriched host and 115 ns for the B-10-enriched host. That calculation retains central-spin dynamics and uses full ab initio hyperfine tensors; its site-resolved inputs cannot be reconstructed from the few longitudinal constants printed in the article.
