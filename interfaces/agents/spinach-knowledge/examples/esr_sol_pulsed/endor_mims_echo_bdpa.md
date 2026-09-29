# examples/esr_sol_pulsed/endor_mims_echo_bdpa.m

Call `endor_mims_echo_bdpa()` with no arguments. This is the stimulated-echo diagnostic stage of the BDPA Mims ENDOR sequence; the nuclear pulse is not applied. The source says the echo becomes sharper as g-tensor anisotropy is increased and estimates a runtime of seconds.

## Spin system and sequence

- The system is one electron and two 1H nuclei at 3.35 T. The diagonal electron g values are 2.00263, 2.00260 and 2.00257. The diagonal electron–proton coupling matrices are `[7.70, 5.30, 2.00] × 10^6` and `[1.00, 1.00, 1.26] × 10^6`; the source does not annotate their units. The basis is `sphten-liouv` without approximation.
- `powder` calls `@endor_mims_echo` in the `esr` context with `rep_2ang_400pts_sph`. τ is 200 ns and the output uses 200 time steps over 0 to 2τ (201 plotted points). Although `n_dur` is assigned 50 μs, this diagnostic stage omits the nuclear pulse; it is not a simulated ENDOR nuclear-pulse response.
- The figure plots `real(answer)` against the constructed time axis in ns and labels intensity in arbitrary units. No data file is saved.

## Requirements and scope

The example needs MATLAB and Spinach, including `powder`, system/basis construction and the `endor_mims_echo` sequence implementation; plotting uses `kfigure`, `kgrid` and standard plotting functions. No DOI is supplied in the source.

Source: [examples/esr_sol_pulsed/endor_mims_echo_bdpa.m](https://github.com/IlyaKuprov/Spinach/blob/main/examples/esr_sol_pulsed/endor_mims_echo_bdpa.m).
