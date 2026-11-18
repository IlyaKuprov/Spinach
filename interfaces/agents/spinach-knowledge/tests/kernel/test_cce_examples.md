# tests/kernel/test_cce_examples.m

[MATLAB source](https://github.com/IlyaKuprov/Spinach/blob/main/tests/kernel/test_cce_examples.m)

## Purpose

Run the CCE application examples and check their independent small-system references.

## Behaviour

- Executes the CCE, gCCE, ME-CCE, and partitioned-CCE examples with figures hidden, then restores the original figure-visibility setting.
- Checks finite numerical results and validity returned by the examples against independently evaluated references in this test.
- Uses observed-spin density matrices and selected-transition signals with the same contracts as the CCE kernel.
- The reduced application models illustrate algorithms; they are not material fits or converged powder calculations.

## Calling convention

```matlab
result=test_cce_examples()
```

**Inputs**

- None; run with Spinach and its test libraries on the MATLAB path.

**Outputs**

- `result` — standard test-harness record containing checks, failure messages, and diagnostics.

## References

- [Registered tests](https://github.com/IlyaKuprov/Spinach/blob/main/tests/lib/test_manifest.m)
- [CCE application examples](../../examples/cce/)

The DEER example check also excludes the pump resonance and requires an empty selected event, guarding frequency selection rather than a fixed physical target.

Powder checks require response/value equality, the averaged final-order response, weighted numerical bounds, and diagnostics consistent with those averaged bounds.

The generalised high-spin examples gcce_sic_clock, gcce_highspin_siv, and gcce_votpp return numeric density arrays and separate clexpn diagnostics. Their independent dense references are local functions of this registered test. Its optional second output, reference_capture, retains the independently constructed physical arrays for numerical-regression capture; the ordinary test-runner first output is unchanged.

The Er:CaWO4, hBN quadrupole, coherent methyl exchange, SiC divacancy, and Si:P donor examples return numeric clexpn answers and separate diagnostics. Their unchanged independent dense bath Hamiltonians, branch propagators, complex amplitudes, and diagnostic arrays are retained in local reference companions and the optional reference_capture output.

The single-nucleus, NV, P1, and partitioned-electron-layer examples return `[answer,diagnostics]` through clexpn and hahn_echo. Analytic/source-pinned/dense references, including the partitioned versus full-order comparison, live in test_cce_examples. The retired entry points are not used.

The native polarised-bath, two-observed-centre, and dissipative-bath examples return [answer,diagnostics]. Their exhaustive sampling comparison, independent dense propagations, negativity, and physical reference arrays live here. The optional reference_capture records those quantities. Assertions explicitly verify the weighted-summation arithmetic bound and transport of positive stock-Liouville arithmetic estimates. Every application example uses the native interface.

The noise example uses native operator descriptors. The cosine reference and original transverse operator are retained in the example-test capture. Native noise tests compare descriptors and explicit Hermitian matrices before connected ensemble reduction and reject invalid descriptor fields.

gcce_finite_drive returns native answer/diagnostics. Its independent two-spin expm propagation and frozen physical reference arrays now live in the test, not in the example.

The migrated cce_deer_bandwidth example returns [answer,diagnostics]. Its independent test uses direct three-spin density-matrix expm propagation with an observed refocusing pulse and selected bath pulse. A native clexpn check moves the spectral window away, verifies the no-pump limit, and verifies allow_empty=false rejects the empty selector. No legacy parameter assembly is used for this example.

## Native molecular powder

cce_molecular_powder returns numeric answer and separate native diagnostics. Its independent full-bath reference assembles rotated physical tensors with Cartesian spin matrices, the stock ordered-pair tensor convention and expm, then averages completed Hahn echoes. Tests cover agreement, the analytic unit initial signal, the final expansion order, and orientation-resolved numerical bounds. There is no legacy entry-point call in any examples/cce file.

Finite-control invocations use finite_drive with native physical descriptors, not a raw sequence field. Independent references retain their physical inputs and assertions; unused raw-event declarations are not part of the reference construction.

Projected Hahn signals use the fixed physical receiver, matching generalised density detection without an extra conjugation of the returned signal. Historical manifold-overlap reference expressions are converted at the reference boundary; pinned external source arrays are retained unchanged.
