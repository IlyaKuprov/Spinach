# examples/nmr_solids/case_studies/cp_square_vs_ramp/cp_square_vs_ramp.m

Source: [cp_square_vs_ramp.m](https://github.com/IlyaKuprov/Spinach/blob/main/examples/nmr_solids/case_studies/cp_square_vs_ramp/cp_square_vs_ramp.m)

## Model and comparison

This example simulates static-powder ¹H-to-¹⁵N cross-polarisation in the doubly rotating frame and compares fixed-amplitude, linear-ramp, and tangent-ramp irradiation. The source describes the comparison as demonstrating the advantages of ramped cross-polarisation and cites [DOI 10.1016/0009-2614(94)00470-6](https://doi.org/10.1016/0009-2614(94)00470-6). The calculation is a simulation; the script supplies no experimental signal for comparison.

The model contains ¹⁵N and ¹H, with sys.magnet set to 9.394, scalar shifts set to zero, coordinates (0, 0, 0) and (0, 0, 1.05), and temperature set to 298. The source does not state units for these assignments. It uses the sphten-liouv formalism with no basis approximation, a 6,400-point spherical powder grid, and the aniso_eq requirement.

## CP sequence and observables

The initial excitation operator is ¹H Lx. The two irradiation operators are ¹H Ly and ¹⁵N Lx, and the receiver coil is ¹⁵N Lx: the reported observable is therefore the ¹⁵N transverse-x expectation following ¹H-to-¹⁵N CP. Each comparison uses 500 time steps of 2 × 10⁻⁶ seconds, giving a 1 ms CP block; the plotted acquisition time axis is labelled in seconds.

The fixed profile assigns 5 × 10⁴ to both irradiation channels for all 500 steps. The linear profile ramps one channel down and the other up over 500 samples. For the tangent profile, the code forms tan(linspace(−1.4, 1.4, 500)), shifts it to start at zero, normalises its maximum to one, then applies the reversed and forward profiles to the two channels at 5 × 10⁴. The source gives no physical unit for these irradiation-power values.

## Output

Spinach powder simulations produce three signals. The figure plots their real time-domain signals and the corresponding ¹⁵N Sx expectation values, with separate traces for constant, linear-ramp, and tangent-ramp CP. The source does not supply measured data or a numerical comparison metric.
