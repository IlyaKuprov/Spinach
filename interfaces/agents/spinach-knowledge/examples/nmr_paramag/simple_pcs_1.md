# examples/nmr_paramag/simple_pcs_1.m

- Signature: `simple_pcs_1()`
- Source: [examples/nmr_paramag/simple_pcs_1.m](https://github.com/IlyaKuprov/Spinach/blob/main/examples/nmr_paramag/simple_pcs_1.m)

## Purpose

Sets up a proton in the presence of a point magnetic-susceptibility centre and computes Redfield relaxation rates. Although the source header mentions pseudocontact shift and Curie relaxation, the executable body reports R1 and R2; it contains no separate PCS calculation or spectrum simulation. The source estimates seconds of calculation time.

## Physical / mathematical content

- Uses a 14.1 T field and a single `1H` spin with scalar shift 2.0 at coordinate [0.0, 0.0, 0.0].
- Places the susceptibility centre at [10.0, 2.5, 3.9] and supplies the symmetric 3x3 susceptibility tensor [[0.0883, -0.0904, 0.0822], [-0.0904, -0.1011, -0.0149], [0.0822, -0.0149, 0.0128]]. The source does not specify units for these tensor and coordinate values.
- Selects Redfield relaxation, zero equilibrium, lab-frame relaxation terms, and a correlation time of `10e-12` s (10 ps).

## Numerical / algorithmic content

The calculation builds a `sphten-liouv` basis with approximation `none`, obtains the relaxation superoperator, and projects it onto the proton `Lz` and `L+` states to calculate R1 and R2. Those two rates are printed; no measured spectrum or experimental comparison is read or produced.

## Implementation structure

- Specifies the proton, field, susceptibility tensor and centre, plus the Redfield parameters.
- Creates the spin system and basis, calls `relaxation`, evaluates the longitudinal and transverse rates from `Lz` and `L+`, and displays them.
