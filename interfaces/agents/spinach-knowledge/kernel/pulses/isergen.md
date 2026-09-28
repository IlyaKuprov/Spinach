# kernel/pulses/isergen.m

- Signature: `H=isergen(HL,HM,HR,dt)`

## Purpose

Second- and fourth-order Iserles product quadrature generators for one time propagation step for a state-independent Hamiltonian.

## Numerical / algorithmic content

- If `HM` is empty, second-order product quadrature is used; otherwise, fourth-order product quadrature is used.

## Parameters / inputs

- `HL` — Hamiltonian at the left edge of the interval.
- `HM` — optional Hamiltonian at the interval midpoint; if empty, second-order quadrature is used.
- `HR` — Hamiltonian at the right edge of the interval.
- `dt` — interval duration, in seconds.

## Outputs

- `H` — effective evolution generator, to be used as `exp(-1i*H*dt)`.
