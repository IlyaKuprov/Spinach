# examples/esr_sol_pulsed/ridme_cu_nitroxide.m

- Signature: `ridme_cu_nitroxide()`

## Purpose

Illustrates RIDME on a Cu(II)–nitroxide two-electron model at `sys.magnet=1.249` (the source describes this as Q-band). The numerical path uses brute-force time propagation and Liouville-space powder averaging, including g-factor orientation effects on the dipolar coupling; extended T1/T2 relaxation is included. The source notes that its analytical calculation uses only isotropic parts of the electron g factors. The source comments give a calculation time of seconds.

## Spin system and relaxation

The two electron spins are ordered Cu first and nitroxide second. Their Zeeman principal values are `[2.056, 2.056, 2.205]` and `[2.009, 2.006, 2.003]`, respectively; both Euler-angle triples are zero. The coordinates are `[43, 0, 0]` and `[0, 0, 0]` Å. Both spins use the `t1_t2` relaxation model: the two R1 entries are `1/(35e-6)` and `1/(2e-3)` Hz, and the R2 entries are `1/(1.5e-6)` and `1/(1.3e-6)` Hz. The model keeps relaxation in the lab frame and sets equilibrium to zero.

The calculation uses the `sphten-liouv` formalism with no basis approximation and disables trajectory-level SSR algorithms. These are choices of this example, not claims about the only valid RIDME setup.

## Sequence, calculation, and output

The initial state is electron `Lz`; spin 2 is the probe. The RIDME timing uses a 16 ns step, 25 and 188 steps in the two evolution dimensions, and a 35 μs mixing time. Powder averaging uses `rep_2ang_800pts_sph`; the ESR calculation is dispatched as `powder(spin_system,@ridme,parameters,'esr')`.

The script forms a phase-cycle trace by summing the real and imaginary components from the PxPxPx, PyPyPx, MxMxPx, and MyMyPx pathways. It then plots the real and imaginary pathway traces together with the corresponding summed RIDME trace against time in μs. The plotted interval is constructed from the two evolution-step counts, so it is not a separately supplied experimental time axis.

## Scope

This is a parameterised simulation example, not a reported experimental fit: the file gives no reference dataset, fitted distance distribution, or expected trace values. Its numerical powder calculation and its stated isotropic-g approximation for the analytical calculation should not be conflated.

Source code: [`examples/esr_sol_pulsed/ridme_cu_nitroxide.m`](https://github.com/IlyaKuprov/Spinach/blob/main/examples/esr_sol_pulsed/ridme_cu_nitroxide.m).
