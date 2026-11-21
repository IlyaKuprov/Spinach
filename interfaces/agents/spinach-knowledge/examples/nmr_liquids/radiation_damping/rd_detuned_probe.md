# rd_detuned_probe: radiation damping

[MATLAB example](../../../../../../examples/nmr_liquids/radiation_damping/rd_detuned_probe.m)

Demonstrates frequency pulling and damping by a phase-shifted probe response.
A small-angle proton signal with a fixed frequency offset and transverse
relaxation is propagated for several feedback phases. The resulting FIDs and
magnitude spectra show how the dissipative and reactive feedback components
change with probe phase.

The small initial tilt keeps the magnetisation near equilibrium. Spectra use
Spinach's first-point correction and zero filling; frequency is in Hz, whereas
operators in the generator use angular frequencies.
