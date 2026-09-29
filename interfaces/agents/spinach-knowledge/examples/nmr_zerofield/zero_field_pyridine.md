# examples/nmr_zerofield/zero_field_pyridine.m

- Signature: `zero_field_pyridine()`
- Source: [`examples/nmr_zerofield/zero_field_pyridine.m`](https://github.com/IlyaKuprov/Spinach/blob/main/examples/nmr_zerofield/zero_field_pyridine.m)

## Purpose

This simulated zero-field NMR example treats 15N pyridine and is set up to reproduce Figure 3 of [Journal of the American Chemical Society, DOI 10.1021/ja2112405](https://doi.org/10.1021/ja2112405). That is the source's stated target, not a report of a run or an independent reproduction check.

## Spin and coupling model

The spin list is five 1H nuclei followed by one 15N nucleus, with the magnetic field set to zero. The source assigns these scalar couplings (spin indices are their positions in that list; values are in Hz): 1–2 and 4–5, 4.88; 1–4 and 2–5, 0.97; 1–3 and 3–5, 1.83; 1–5, −0.12; 2–3 and 3–4, 7.62; 2–4, 1.38; 1–6 and 5–6, −10.93; 2–6 and 4–6, −1.47; 3–6, 0.27; and 6–6, 0.00. The Hilbert-space basis is `zeeman-hilb` with no approximation. The example does not specify spatial coordinates, gradients, chirps, or a field-drop schedule.

## Acquisition and processing

The FID is calculated by `liquid(...,@zerofield,...,'labframe')`. The settings give a 60 Hz sweep, 512 acquired points, 1024-point zero filling, zero offset, 1H excitation, a π/2 flip angle, uniaxial detection, and an axis in Hz. The code subtracts the FID mean, applies exponential apodisation with parameter 12, Fourier transforms and shifts the signal, then plots the real spectrum. The source comment estimates calculation time as seconds; that is a source note, not a runtime measured here.

This is a zero-field NMR simulation, not SPEN, ultrafast DOSY, or a multiple-quantum experiment: the source has no spatial encoding or gradient schedule and no explicit multiple-quantum selection. It imports no experimental data.
