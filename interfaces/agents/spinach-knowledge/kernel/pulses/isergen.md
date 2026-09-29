# kernel/pulses/isergen.m

[Source on GitHub](https://github.com/IlyaKuprov/Spinach/blob/main/kernel/pulses/isergen.m) · [Spin Dynamics Wiki: isergen.m](https://spindynamics.org/wiki/index.php?title=isergen.m)

Signature: `H=isergen(HL,HM,HR,dt)`

## Purpose and interval data

Builds the effective generator for one time-propagation interval with a state-independent Hamiltonian. `HL` and `HR` are the Hamiltonians at the left and right edges; `HM` is optional and, when nonempty, is the midpoint Hamiltonian. `dt` is the interval duration in seconds. The output is used in `exp(-1i*H*dt)`.

## Quadrature choice and ordering

An empty `HM` selects the second-order product quadrature:

`H=(HL+HR)/2 + (1i*dt/6)*(HL*HR-HR*HL)`

A supplied `HM` selects the fourth-order product quadrature:

`H=(HL+4*HM+HR)/6 + (1i*dt/12)*(HL*HR-HR*HL)`

The commutator correction is ordered as `HL*HR-HR*HL`; reversing its factors changes the expression. These are endpoint samples, plus a midpoint sample for the fourth-order option—not a waveform file or a sequence of pulse-amplitude samples. Any pulse amplitude or phase is represented through the supplied Hamiltonians; there are no separate pulse-control arguments. Plotting and file output are not handled here.

The source checks that each supplied Hamiltonian is square and that `dt` is a real numeric scalar. It does not explicitly check that the Hamiltonian dimensions match one another.
