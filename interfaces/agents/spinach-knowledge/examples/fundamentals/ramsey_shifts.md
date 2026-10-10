# examples/fundamentals/ramsey_shifts.m

- MATLAB implementation: [examples/fundamentals/ramsey_shifts.m](https://github.com/IlyaKuprov/Spinach/blob/main/examples/fundamentals/ramsey_shifts.m)

- Signature: `ramsey_shifts()`
- Source: `examples/fundamentals/ramsey_shifts.m`

## Purpose

Shows how an off-resonant proton-channel drive shifts the phases of 13C and 15N in a three-spin system. For nucleus n, the signed Ramsey frequency shift is

`delta_n = w_n w1n^2 / (w_n^2 - w_c^2)`,

where w_n is the signed Zeeman frequency, w_c the signed carrier frequency, and w1n the proton-drive amplitude scaled by the nucleus-to-proton magnetogyric-ratio ratio. Because 15N has a negative magnetogyric ratio, its phase evolves in the opposite direction to 13C.

## Model and comparison

The script forms the transverse states of 13C and 15N in a 1H/13C/15N system and propagates their normalised sum for 20 ms under a constant proton drive. The off-resonant nuclei are not driven directly. A Ramsey/Bloch–Siegert correction is enabled for the proton-channel control, and the final accumulated phases are compared with the analytic shifts using the signed base frequencies and magnetogyric ratios.

At 14.1 T and a drive amplitude of 2π·25 kHz, the example checks each numerical phase against its analytic value with relative tolerance 10^-6 and checks that the 15N and 13C phases have opposite signs. It then doubles the amplitude and requires the 13C phase ratio to be 4 (quadratic drive-amplitude scaling), and halves the field and requires a ratio of 2 (inverse-field scaling), each within 10^-5.

The zero drift uses the compiled state-space dimension `bas.offsets(end)`, matching the state and control-operator dimensions.
