# examples/esr_sol_pulsed/hard_3_pulse_deer_cu.m

- Signature: `hard_3_pulse_deer_cu()`
- Source: [`hard_3_pulse_deer_cu.m`](https://github.com/IlyaKuprov/Spinach/blob/main/examples/esr_sol_pulsed/hard_3_pulse_deer_cu.m)

## Purpose and spin model

This example compares a three-pulse DEER time trace for a Cu(II)–NO pair with an analytical trace. The numerical calculation uses brute-force time propagation and powder averaging in Liouville space, retaining the orientation dependence of the electron g tensors in the dipolar interaction. The analytical comparison replaces each g tensor by its isotropic mean. The source describes the setup as X-band and lists a static field of `0.33` T.

There are two electron spins (`E`): spin 1 has g-tensor principal values `[2.056, 2.056, 2.205]`, and spin 2 has `[2.009, 2.006, 2.003]`; both Euler-angle triples are `[0, 0, 0]`. Their coordinates are `[0, 0, 0]` and `[20, 0, 0]` Å, so the pair is separated along x by 20 Å. The pairwise dipolar coupling is obtained from those coordinates by `xyz2dd`; no exchange term is set. Spinach uses the untruncated `sphten-liouv` basis (`bas.approximation='none'`) and disables trajectory-level SSR.

## Pulse sequence and sampling

The three-pulse sequence is delegated to `deer_3p_hard_deer` (the [shared sequence helper](https://github.com/IlyaKuprov/Spinach/blob/main/experiments/esr_dipolar/deer_3p_hard_deer.m)), called through `powder` in the `'deer'` context with brief output. The helper applies a hard `π/2` probe pulse, evolves the trajectory for the configured interval, applies a hard `π` pump pulse, refocuses the stored trajectory, applies a hard `π` probe pulse, then records final evolution through the probe detection channel. Here the interval is 50 × 10 ns = 500 ns; the time axis has 51 samples from 0 through 500 ns. Powder averaging uses `rep_2ang_1600pts_sph`. The caller sets no finite pulse widths or independent offsets.

## Calculated traces and figure

The numerical observable is `deer.deer_trace`. For the comparison, the source calculates the dipolar coupling with `xyz2dd`, scales it by the product of the two mean g values divided by `spin_system.tols.freeg^2`, then evaluates `0.35*deer_analyt(D,0,time_axis)`. A figure places `imag(deer_num.deer_trace)` beside that analytical result; both are plotted against time in microseconds, with the displayed analytical axes limited to 0–0.5 μs and −0.1–0.4. The source estimates the calculation time in seconds. The script creates a figure but does not save a data file or figure.

## Source

[Spinach example source](https://github.com/IlyaKuprov/Spinach/blob/main/examples/esr_sol_pulsed/hard_3_pulse_deer_cu.m)
