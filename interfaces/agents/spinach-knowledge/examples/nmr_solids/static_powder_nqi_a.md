# examples/nmr_solids/static_powder_nqi_a.m

- Signature: `static_powder_nqi_a()`
- Source: [examples/nmr_solids/static_powder_nqi_a.m](https://github.com/IlyaKuprov/Spinach/blob/main/examples/nmr_solids/static_powder_nqi_a.m)

## Purpose and model

Calculates the static powder 14N NMR pattern of L-valyl-L-alanine. The source says the large orientation grid is intended to reproduce Figure 5 of O'Dell and Ratcliffe ([DOI](https://doi.org/10.1016/j.cplett.2011.08.030)); its runtime estimate is minutes.

Two 14N spins are assigned separate quadrupolar interactions with `eeqq2nqi`: values 1.24e6 and 3.06e6, asymmetries 0.22 and 0.40, respectively. The field parameter is 21.1. A full Zeeman Hilbert-space basis is used with no approximation. `powder` performs the static orientation average on `icos_2ang_163842pts`; no rotor or gradient parameters are set.

## Acquisition and processing

The NMR acquisition selects 14N with offset 0, sweep 6e6, 512 acquired points, and 2048-point zero-fill. The frequency-axis setting is MHz and the axis is inverted. The initial and detection states are both 14N `L+`. The powder FID is exponentially apodised with parameter 6, Fourier transformed, and plotted using its real part.
