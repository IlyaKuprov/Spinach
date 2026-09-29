# examples/esr_sol_pulsed/hard_3_pulse_deer_no.m

- Signature: `hard_3_pulse_deer_no()`

## Purpose

Models a three-pulse double electron–electron resonance (DEER) experiment on a pair of nitroxide radicals at X-band. The two electron spins are placed 25 Å apart. The calculation uses brute-force time propagation and powder averaging, so anisotropic g tensors and the orientation-dependent dipolar interaction contribute to the orientational average.

## Spin system and fixed inputs

The field is 0.33 T; both spins are electrons (E) with principal g values [2.0089, 2.0061, 2.0027]. Their Euler-angle inputs are [1, 2, 3] and [3, 1, 2]; the example does not label the angle units. The coordinates are [0, 0, 0] and [25, 0, 0] Å. The calculation uses the exact Zeeman-Hilbert basis (zeeman-hilb, approximation none). These are one fixed model, not a sweep.

## Pulse sequence and sampled signal

The script prepares Lz magnetisation and detects the probe-spin response using the spin-1 coil state. Selective transverse Lx operators address spins 1 and 2 as probe and pump. It calls the three-pulse hard-DEER helper through powder averaging. The helper uses ideal selective pulses: probe π/2, a probe evolution trajectory, pump π with refocusing evolution, probe π, and final evolution. Here the trajectory has 100 increments of 10 ns, for a 1 µs time axis with 101 plotted points. The brief-output setting requests the DEER trace rather than the helper's optional detailed pulse FIDs.

The powder average uses the rep_2ang_3200pts_sph spherical grid. The plotted quantity is the imaginary part of deer_trace versus time in microseconds. The script creates a figure but specifies no saved data file; its comment estimates seconds for calculation time.

## Reference and implementation

The example gives the nitroxide g-tensor reference as [DOI: 10.1063/1.1697233](http://dx.doi.org/10.1063/1.1697233). See the [example source](https://github.com/IlyaKuprov/Spinach/blob/main/examples/esr_sol_pulsed/hard_3_pulse_deer_no.m) and the [three-pulse DEER helper](https://github.com/IlyaKuprov/Spinach/blob/main/experiments/esr_dipolar/deer_3p_hard_deer.m).
