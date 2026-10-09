# experiments/nmr_protein/hnco.m

[Canonical MATLAB source](https://github.com/IlyaKuprov/Spinach/blob/main/experiments/nmr_protein/hnco.m) · [Spinach Wiki](https://spindynamics.org/wiki/index.php?title=hnco.m)

## What the sequence represents

This is a phase-sensitive, three-dimensional HNCO sequence for 1H, 13C, and 15N-labelled proteins. Its frequency dimensions are F1 = N, F2 = CO, and F3 = H. The source cites the reported HNCO experiment ([DOI: 10.1016/0022-2364(90)90333-5](https://doi.org/10.1016/0022-2364(90)90333-5)) and the bidirectional-propagation method ([DOI: 10.1016/j.jmr.2014.04.002](https://doi.org/10.1016/j.jmr.2014.04.002)).

The implementation starts from concentration-weighted NH-proton longitudinal magnetisation through `state` with the `cheap` method, and uses unweighted NH-proton `L+` from `coil_state(spin_system,'L+',find(HNs),'cheap')` as the detection state. Spin selections use PDB atom labels: the code selects `H`, `N`, `C` (the carbonyl carbon), and `CA` from `spin_system.comp.labels`; isotope-wide pulse operators are requested as `1H` and `15N`. The source comment requires PDB atom IDs such as `CA` and `HA` to be represented in `sys.labels`. The pulse operators are built from `L+` and converted to Cartesian x/y operators for ideal `step` pulses.

The forward/reverse trajectories implement the H-to-N-to-CO correlation pathway and the three acquisition dimensions. The source separately selects positive and negative N coherence for F1 and positive/negative carbonyl coherence for F2; `stitch` combines those forward and backward pathways for States quadrature. `parameters.f1_decouple` selects the source's T1 refocusing branch: its enabled branch inverts H together with CO and CA, while the other branch inverts CO and CA. This is the implemented switch, not a general-purpose decoupling waveform.

## Inputs and timing

Call signature: `fid=hnco(spin_system,parameters,H,R,K)`.

- `parameters.npoints`: three integers ordered `[t1 t2 t3]`.
- `parameters.sweep`: three sweep widths in Hz, ordered `[f1 f2 f3]`.
- `parameters.tau`: three delays in seconds. The source gives `[2.25e-3, 14e-3, 4e-3]` as reasonable values.
- `parameters.f1_decouple`: logical 0/1 switch described above.
- `H`, `R`, `K`: Hamiltonian, relaxation, and kinetics matrices supplied by the context function.

The implementation requires the `sphten-liouv` formalism and matrix-valued `H`, `R`, and `K` with matching dimensions. It checks for three entries in each sampling/timing vector and accepts `f1_decouple` as 0 or 1. Initial and detection states are created inside the routine; this function does not expose `rho0` or `coil` parameters.

## Output

`fid` has four fields: `pos_pos`, `pos_neg`, `neg_pos`, and `neg_neg`, the paired coherence-sign components for States quadrature. Each stitched 3D array is permuted with `[3 2 1]`, so its dimension order is `[t3,t2,t1]` (nominal extents `[npoints(3),npoints(2),npoints(1)]`); the corresponding frequency labels are H, CO, and N.
