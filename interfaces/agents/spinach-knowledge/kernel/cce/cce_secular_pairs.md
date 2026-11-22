# kernel/cce/cce_secular_pairs.m

[`cce_secular_pairs`](../../../../../kernel/cce/cce_secular_pairs.m) computes a secular CCE-2 Hahn echo for an unpolarised bath of identical spin-½ nuclei with longitudinal dipolar hyperfine couplings. This direct helper and the compact `clexpn` Hahn route share the streamed pair reduction. The Er:CaWO₄ demonstration uses the native context with a compact created system.

```matlab
echo=cce_secular_pairs(xyz,field_axis,hyperfine,dipolar,times,abs_tol)
```

Supply distinct centre-relative sites in metres as an N-by-3 array, a 3-by-1 field unit vector, the two dipolar prefactors in rad s⁻¹ m³, and a column of total echo times in seconds. The returned column contains dimensionless normalised amplitudes. Each single-spin contribution is unity. The two-spin contribution comes from the secular flip-flop block; central and common nuclear Zeeman phases refocus.

`abs_tol` bounds the total loss omitted by pair screening: each rejected pair has a loss upper bound no greater than its share of the budget. Setting `abs_tol=0` omits only pairs whose loss bound is zero. The registered `kernel/cce_erbium_s12` test compares the pair product with dense Hahn propagation and checks unity for a singleton and for vanishing hyperfine detuning.
