# kernel/optimcon/bss_ops.m

Source: [kernel/optimcon/bss_ops.m](https://github.com/IlyaKuprov/Spinach/blob/main/kernel/optimcon/bss_ops.m)

## Purpose

Builds one Bloch–Siegert response operator for each control channel. In GRAPE, the coefficient multiplying a channel's response operator in a time slice is the square of that channel's physical control amplitude. The routine constructs the operators; it does not itself calculate an objective or an objective gradient.

For spin n and channel c, define the signed laboratory frequency <code>omega_n = spin_system.inter.basefrqs(n)</code> and signed carrier <code>omega_c = carrier_frq(c)</code>. The coefficient multiplying that spin's longitudinal operator is <code>(gamma_n/gamma_c)^2/[2*(omega_n+omega_c)]</code> for every spin, plus <code>(gamma_n/gamma_c)^2/[2*(omega_n-omega_c)]</code> for foreign-isotope spins. The on-channel resonant term is omitted because it is the control operator that GRAPE propagates exactly. Only spin-type particles marked <code>S</code> are included.

## Call and data

<code>resp_ops = bss_ops(spin_system,channels,carrier_frq)</code>

- <code>spin_system</code> supplies the spin-system formalism, composition, isotope data and signed base frequencies.
- <code>channels</code> is a cell array of character vectors naming the control-channel isotopes. There is one output operator per channel; for an X/Y pair, for example, the channel list may be <code>{'1H','1H'}</code>.
- <code>carrier_frq</code> supplies one signed carrier frequency per channel, in rad/s. The documentation specifies a row vector; validation checks the number of elements, not row orientation.
- <code>resp_ops</code> is a cell array of the response operators in the supplied spin-system formalism.

There are no default input values. The routine requires <code>spin_system.comp</code>; it checks that <code>channels</code> is a cell array of character vectors naming isotopes in <code>spin_system.comp.isotopes</code>. Carrier frequencies must be numeric, real, finite, nonzero, and one per channel. It rejects a spin when <code>abs(omega_n+omega_c) &lt; 1e-6*abs(omega_c)</code>; for foreign isotopes it also rejects <code>abs(omega_n-omega_c) &lt; 1e-6*abs(omega_c)</code>.

The formulas use angular-frequency values in rad/s. No gradient or adjoint is returned.
