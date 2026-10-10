# examples/kinetics/nonlinear/diels_alder_spec.m

Source: https://github.com/IlyaKuprov/Spinach/blob/main/examples/kinetics/nonlinear/diels_alder_spec.m

- Signature: `diels_alder_spec()`

## Model and scope

This script couples a second-order Diels-Alder model, `acetylene (A) + butadiene (B) -> cyclohexadiene (C)`, to proton spin dynamics. It imports A, B, and C from `acetylene.out`, `butadiene.out`, and `cyclohexadiene.out`; each `g2spinach` call passes 31.8 as its third argument (the script does not state its unit), with a minimum imported J coupling of 2.0 Hz. Ethanol (D) is a six-proton solvent subsystem; the script does not assign it coordinates. The source assigns three 7.0 Hz couplings from ethanol protons 1-3 to protons 4-5. The groups occupy spins 1-2, 3-8, 9-16, and 17-22, respectively.

The concentration model starts at `[0.01, 0.02, 0, 17.1] mol/L` for `[A, B, C, D]`. The explicit additive A+B→C record has rate 25.0 L/(mol*s) and the original eight-spin atom matching. True concentrations are stored in `chem.concs`; tracing spins gives the four-pool concentration-only kernel model. D remains a spectator. Concentrations are advanced over 0-10 s in 100 LG4 steps and plotted for A, B, and C only. For the spin model, the field is 14.1 T; Redfield and T1/T2 relaxation are configured with secular retention and zero equilibrium. Correlation times are 1e-12, 20e-12, 50e-12, and 5e-12 s for A-D; solvent R1 and R2 entries are set to 0.5.

## Acquisition and output

The concentration-weighted spin trajectory combines `unit_state` populations and the weighted A–C longitudinal preparation, then uses the kernel reaction generator and a two-point Lie step. Interpolated populations enter through unit coordinates; additive product unit arrival is shared equally between reactants. `coil_state` detects without a second concentration factor. Nine 1H pulse-acquire simulations start at integer times 0-8 s; each applies a pi/2 Ly pulse and evolves 4096 points at 4000 Hz with offset 2370 Hz. The script assembles the sparse history-dependent chemistry on the CPU and transfers each interval-edge generator with the evolution arrays to the GPU, applies exponential apodisation with parameter 6, zero-fills to 16384 points, and plots a waterfall of real spectral intensity versus chemical shift (ppm) and start time (s).

This describes the configured example, not an independently established reaction yield or experimental spectrum. The source header estimates hours of calculation and notes hard-coded GPU use.
