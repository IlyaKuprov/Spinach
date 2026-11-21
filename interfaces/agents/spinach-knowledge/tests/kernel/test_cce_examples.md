# tests/kernel/test_cce_examples.m

[MATLAB source](https://github.com/IlyaKuprov/Spinach/blob/main/tests/kernel/test_cce_examples.m)

## Purpose

Check compact CCE application models against independent small-system references. The test does not run public figure examples.

## Behaviour

- Discovers public MATLAB functions in `examples/cce` and checks that each takes no arguments and returns no results. Private helpers are excluded and the figure calculations are not run.
- Executes CCE, gCCE, ME-CCE, and partitioned-CCE regression fixtures without plotting or needing the public examples on the MATLAB path.
- Checks finite numerical results and validity against independently evaluated references.
- Uses observed-spin density matrices and selected-transition signals with the same contracts as the CCE kernel.
- The reduced application models illustrate algorithms; they are not material fits or converged powder calculations.

## Calling convention

```matlab
[result,reference_capture]=test_cce_examples()
```

**Inputs**

- None; run with Spinach and its test libraries on the MATLAB path.

**Outputs**

- `result` — standard test-harness record containing checks, failure messages, and diagnostics.
- `reference_capture` — independent physical arrays and numerical diagnostics keyed by the original case names. The ordinary test runner needs only `result`.

## References

- [CCE application examples](../../examples/cce/)

The DEER check compares direct three-spin density-matrix propagation with a selected bath pulse and observed-spin refocusing. Moving the spectral window away verifies the no-pump limit and verifies that `allow_empty=false` rejects the empty selector.

Powder checks require response/value equality, the averaged final-order response, weighted numerical bounds, and diagnostics consistent with those averaged bounds.

The generalised SiC clock, high-spin SiV, and VO(TPP) fixtures compare density arrays with independent dense references. The Er:CaWO₄, hBN quadrupole, coherent methyl exchange, SiC divacancy, and Si:P donor fixtures retain their dense bath Hamiltonians, branch propagators, complex amplitudes, and diagnostics for reference capture.

The single-nucleus, NV, P1, and partitioned-electron-layer fixtures exercise `clexpn` and `hahn_echo`, with analytic, pinned external, or dense references. The electron-layer case compares partitioned against full-order propagation.

The polarised-bath fixture compares sampled against exhaustive averaging and checks the weighted-summation arithmetic bound. The two-observed-centre and dissipative-bath fixtures check dense propagation and density physicality, including negativity for the two centres and positive stock-Liouville arithmetic estimates for the dissipative bath.

The noise fixture uses native operator descriptors and captures the analytic cosine reference and original transverse operator. The separate native noise tests compare descriptors and explicit Hermitian matrices before connected ensemble reduction and reject invalid descriptor fields.

The finite-drive fixture uses native physical control descriptors and checks independent two-spin `expm` propagation.

## Native molecular powder

The powder fixture assembles rotated physical tensors using Cartesian spin matrices and the stock ordered-pair tensor convention. Independent full-bath `expm` propagation is performed before averaging the completed Hahn echoes. Tests cover agreement, the analytic unit initial signal, the final expansion order, and orientation-resolved numerical bounds.

Projected Hahn signals use the fixed physical receiver, matching generalised density detection without an extra conjugation of the returned signal. Historical manifold-overlap reference expressions are converted at the reference boundary; pinned external source arrays are retained unchanged.
