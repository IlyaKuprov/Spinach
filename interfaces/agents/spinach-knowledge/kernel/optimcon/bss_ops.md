# kernel/optimcon/bss_ops.m

- Signature: `resp_ops=bss_ops(spin_system,channels,carrier_frq)`

## Purpose

Builds one Bloch–Siegert response operator for each optimal-control channel. Its coefficient in each GRAPE time slice is the square of that channel's physical control amplitude. For spin `n) and channel `c`, the second-order frequency-shift coefficient is

```
(gamma_n/gamma_c)^2 * [1/(2*(omega_n+omega_c)) + 1/(2*(omega_n-omega_c))]
```

The second term is included only for spins whose isotope differs from the channel isotope. Thus an on-channel spin receives only the never-resonant term; a foreign-isotope spin receives both terms. The coefficient multiplies that spin's longitudinal operator. The signed frequencies `omega_n` and `omega_c` are the spin's laboratory-frame Zeeman frequency and the channel carrier frequency, respectively.

## Syntax

```matlab
resp_ops=bss_ops(spin_system,channels,carrier_frq)
```

## Inputs

- `spin_system` — Spinach spin-system structure.
- `channels` — cell array of isotope strings, one per control channel (for example, `{'1H','1H'}` for an X/Y pair).
- `carrier_frq` — row vector of signed carrier frequencies in rad/s, one per channel. For an on-resonance transmitter, use the corresponding `spin_system.inter.basefrqs` value.

## Output

- `resp_ops` — cell array of Bloch–Siegert response operators, one per channel, in the spin-system formalism.

The routine requires each channel isotope to be present in the spin system and rejects zero carriers and degenerate frequency denominators.
