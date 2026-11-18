# examples/cce/gcce_highspin_siv.m

- Signature: `[answer,diagnostics]=gcce_highspin_siv()`

## Purpose

Spin-3/2 vacancy-like centre with mixed spin-half nuclear environment. The independent dense Kronecker/expm comparison is retained in test_cce_examples.

## Spin system and sequence

Standard Spinach sys/inter inputs are passed to create; clexpn compiles metadata without creating a full-bath basis. The example is self-contained, takes no inputs, uses SI times and explicit Hz coupling tensors, and opens a plot. It uses a local one-worker pool. Spinach and EasySpin must not share the MATLAB path.

## Observable and limits

The numeric answer contains all observed-spin density elements in column-major order, one row per time. The separate diagnostics retain validity masks and conservative numerical bounds; no clipping or positivity repair is applied. Independent reference checks are in test_cce_examples.

This is a deliberately small illustrative system, not a publication-scale material reproduction. Material structures, tensor frames, defect origins, and source licences must be resolved before quantitative experimental comparisons. A finite bath does not establish a thermodynamic T2 or microscopic T1. Read the source header in [gcce_highspin_siv.m](https://github.com/IlyaKuprov/Spinach/blob/main/examples/cce/gcce_highspin_siv.m).
