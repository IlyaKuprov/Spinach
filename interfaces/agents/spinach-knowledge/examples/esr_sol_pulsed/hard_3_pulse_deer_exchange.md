# examples/esr_sol_pulsed/hard_3_pulse_deer_exchange.m

- Signature: `hard_3_pulse_deer_exchange()`
- Source: [`hard_3_pulse_deer_exchange.m`](https://github.com/IlyaKuprov/Spinach/blob/main/examples/esr_sol_pulsed/hard_3_pulse_deer_exchange.m)

## Physical aim and model

The example calculates a three-pulse DEER trace for two Cu(II) electron spins in a linked porphyrin complex with strong exchange coupling. It averages the trace over a prescribed distribution of exchange couplings, using brute-force time propagation and numerical powder averaging in Liouville space.

The field is `1.2132` T. Both `E` spins have the same Zeeman tensor, `[2.050, 2.050, 2.195]`, with Euler angles `[0, 0, 0]`. Their coordinates are `[0, 0, 0]` and `[24.50, 0, 0]` Å. The scalar exchange coupling is scanned over 20 values from `6e6` to `20e6` Hz. The code assigns each point a weight with `gaussfun(j_values-13.1e6,4.2622e6)` and normalises the weights to sum to one; it does not label the second `gaussfun` argument as a standard deviation. For each coupling value the spin system is rebuilt with that scalar exchange term and an untruncated `sphten-liouv` basis.

## Three-pulse protocol and numerical observable

The `deer_3p_hard_deer` [sequence helper](https://github.com/IlyaKuprov/Spinach/blob/main/experiments/esr_dipolar/deer_3p_hard_deer.m) applies ideal hard pulses: `π/2` on the probe, trajectory evolution for the configured interval, `π` on the pump, refocusing evolution, then `π` on the probe followed by final evolution and probe-channel detection. It is called through `powder` in the `'deer-zz'` context with brief output and the `rep_2ang_400pts_sph` grid. This context omits flip-flop terms of the inter-electron dipolar interaction. The interval is 200 × 2.5 ns = 500 ns (`stepsize=1e-8/4` s), and the trace has 201 samples over 0–500 ns; no finite pulse widths or separate offsets are configured.

The observable from each calculation is `deer.deer_trace`; the script forms a normalised weighted sum of the traces across the 20 exchange values. It plots `imag(answer)` against time in microseconds. The source estimates the calculation time in minutes. There is no analytical comparator or file-writing operation in the script; its output is a figure.

## Source

[Spinach example source](https://github.com/IlyaKuprov/Spinach/blob/main/examples/esr_sol_pulsed/hard_3_pulse_deer_exchange.m)
