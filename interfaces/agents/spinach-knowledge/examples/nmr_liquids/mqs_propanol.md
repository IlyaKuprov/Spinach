# examples/nmr_liquids/mqs_propanol.m

- Signature: `mqs_propanol()`
- Source: [`examples/nmr_liquids/mqs_propanol.m`](https://github.com/IlyaKuprov/Spinach/blob/main/examples/nmr_liquids/mqs_propanol.m)

## Purpose and model

A simulated multiple-quantum NMR experiment for the source's seven-proton propanol spin system; the source estimates seconds of calculation time. The field parameter is 14.1. The chemical-shift entries are `[3.438, 3.438, 1.429, 1.429, 0.775, 0.775, 0.775]`. The 2J coupling assignments are 7.5 for pairs (1,3), (1,4), (2,3), and (2,4), and 7.0 for (3,5), (3,6), (3,7), (4,5), (4,6), and (4,7). The 3J assignments are 0.5 for each pair of spins 1 or 2 with spins 5, 6, and 7. The source does not annotate units for these entries. The full `sphten-liouv` basis uses no approximation and exploits `S2`, `S2`, and `S3` permutation groups on spin sets `[1 2]`, `[3 4]`, and `[5 6 7]`.

## Coherence experiment and processing

The initial state is proton `Lz` and the detected state is proton `L+`. The sequence parameters set an angle of `pi/2`, offsets `[1200 1200]`, sweeps `[6000 2700]`, 512 points and 2048 zero-filled points in each dimension, two 1H dimensions, and kHz axis units. The offset and sweep units are not stated in the source. The selected coherence order is +3.

The `@mqs` liquid-NMR sequence is simulated at three tau delays: 0.0333, 0.0710, and 0.5000 s. Each 2D FID is sine-apodised in both dimensions, transformed with a zero-filled 2D FFT, and displayed as an absolute-value spectrum. Thus the three panels compare the selected +3 coherence response across tau; the source does not supply experimental intensities or a fitted transfer rate. No DOI is cited.
