# examples/quantum_tech/diamond_defects/diamond_co_epr_xw.m

- Source: [diamond_co_epr_xw.m](https://github.com/IlyaKuprov/Spinach/blob/main/examples/quantum_tech/diamond_defects/diamond_co_epr_xw.m)
- Signature: `diamond_co_epr_xw()`

## Purpose and model

Simulate field-swept powder EPR spectra for the O4 cobalt-related centre in diamond at X and W bands. The example passes `centre='o4'` and `orientation='111'` to `diamond_co`, then selects the electron signal with `parameters.spins={'E'}`. The builder supplies the centre's magnetic tensors; the example does not define a new defect Hamiltonian. Its cited magnetic-parameter source is Nadolinny et al., *Crystals* **7**, 237 (2017), [doi:10.3390/cryst7080237](https://doi.org/10.3390/cryst7080237).

The calculation uses the `zeeman-hilb` basis with no basis approximation and the source's powder grid `rep_2ang_100pts_sph`. It sets `sys.magnet=1`, `fwhm=1e-4`, `int_tol=1e-5`, `tm_tol=0.1`, `npoints=1024`, and `rspt_order=Inf`; these are simulation inputs, not experimentally fitted values in this example.

## Sweep sequence and output

The script runs two independent field sweeps on the same spin system:

| Band | Microwave frequency | Field window |
| --- | ---: | ---: |
| X | 9.5 GHz (`9.5e9`) | 0.2–0.45 T |
| W | 94 GHz (`94e9`) | 2.6–4.2 T |

The returned `spec_x` and `spec_w` spectra are plotted against their respective `b_axis` values in tesla; the plotted intensity is labelled in arbitrary units. These are simulated spectra from the spin model, not new EPR measurements. The source comment estimates minutes for the calculation; that estimate was not benchmarked here.
