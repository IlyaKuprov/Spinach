# examples/nmr_liquids/inad_three_spin_2d.m

- Signature: `inad_three_spin_2d()`
- Source: [examples/nmr_liquids/inad_three_spin_2d.m](https://github.com/IlyaKuprov/Spinach/blob/main/examples/nmr_liquids/inad_three_spin_2d.m)

## Purpose

A two-dimensional liquid-state INADEQUATE example for a generic three-spin `13C` system. The indirect F1 axis is the double-quantum dimension and F2 is detected in the one-quantum dimension. The source estimates calculation time in seconds and credits Theresa Hune and Christian Griesinger.

## Implementation

At `16.44` T (700 MHz), the model has three `13C` spins with shifts `10`, `30`, and `70` ppm; its nonzero couplings are `J(1,2)=20` Hz and `J(1,3)=60` Hz, with no coupling between spins 2 and 3. It uses the full `sphten-liouv` basis with no approximation and enumerates pair-labelled `13C` isotopomers.

The two-dimensional sequence observes `13C` without decoupling, uses `J=50` Hz and offset value `17604.78` (the source comment places the offset at 100 ppm), and specifies sweep values `[35213.086 34722.223]` with a source comment identifying the sweep width as 200 ppm. The acquisition grid is `[128 2048]` points, zero-filled to `[512 8192]`; the source does not assign units to the numeric offset or sweep-vector entries. For each isotopomer, `inadequate_2d` returns cosine and sine FIDs. Both receive cosine apodisation in both dimensions; F2 is Fourier transformed separately, the real parts are recombined as a States signal, and F1 is Fourier transformed. The plotted real spectrum labels F1 as the DQ dimension in ppm.
