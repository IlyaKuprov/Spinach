# tests/kernel/test_nmr_liquids_alignment_suite.m

- Signature: `result=test_nmr_liquids_alignment_suite()`

## Purpose

Regression test for compact liquid-state NMR pulse-sequence paths. Returns a test result with explanatory messages.

## Checks

- On a two-`1H` system, runs gCOSY `P`, `N`, and `P+N` gradient pathways. Checks that the separate FIDs are finite, non-zero, and distinct; that `P+N` returns `pos` and `neg` branches matching the separate acquisitions; and that Fourier recombination produces a finite spectrum with non-zero real signal.
- On a coupled `1H`–`15N` system, checks that HMBC with `15N` as the selected heteronucleus produces a finite, non-zero FID. Checks finite, non-zero positive and negative pathways for HSQC and CT-HSQC, including the expected CT-HSQC array sizes.
- Runs NOESY-HSQC through its decoupled mixing path and checks that all four pathway components are finite.
- On a single-`1H` system, compares direct TOCSY outputs at zero and non-zero mixing times and checks that relaxation attenuates the combined `cos` and `sin` signal during the spin-lock interval.