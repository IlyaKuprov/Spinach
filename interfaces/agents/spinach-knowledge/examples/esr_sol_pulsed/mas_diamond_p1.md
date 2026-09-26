# examples/esr_sol_pulsed/mas_diamond_p1.m

- Signature: `mas_diamond_p1()`

## Purpose

Two-pulse echo-detected frequency-swept EPR spectra of the P1 substitutional nitrogen defect in diamond, static and under magic angle spinning, after Figure 1a of Khamrui et al., J. Phys. Chem. Lett. 2026, <https://doi.org/10.1021/acs.jpclett.6c02108>: 400 ns long pulses with 416 kHz nutation frequency, 300 ns interpulse delay, at 6.9 T, static and at 10, 25, and 37 kHz MAS. The carrier is stepped across the spectrum and the integrated echo is recorded at each frequency, with the spectra normalised to the static one. Calculation time: hours on a 256-core node.

## Physical / mathematical content

- Pulsed ESR / EPR solid-state examples. These scripts revolve around electron spin echo sequences, DEER, RIDME, ENDOR, ESEEM, and HYSCORE. They combine anisotropic Zeeman and hyperfine Hamiltonians with selective pulses, echo formation, and orientation averaging.
- The 14N hyperfine coupling of the P1 centre (dipolar part 10.9 MHz) makes the two outer lines dephase under spinning as their resonance frequencies move during the sequence, whereas the central line survives.
- All P1 centre tensors are axial and coaxial, so a two-angle powder grid is sufficient; relaxation is omitted: with T2 long against the echo window it scales the four spectra by nearly the same factor, which the normalisation removes.

## Known limitations against Figure 1a

- Line positions and the collapse of the outer lines under spinning match the paper, but the intensities do not reproduce the paper's own simulation: normalised to the static central peak, the static outer perpendicular edges come out at 0.20 and 0.18 (paper about 0.4), and the central line keeps 0.98, 0.90, and 0.82 of its echo at 10, 25, and 37 kHz (paper about 0.87, 0.52, and 0.33). The paper's simulated outer lines are broad humps where this example gives the perpendicular-edge singularity convolved with the roughly 1 MHz pulse response.
- Zeroing the 14N quadrupole leaves the 37 kHz central survival unchanged, so it is not the source. Untested candidates: the paper's model keeps only the secular hyperfine term with no nuclear Zeeman interaction, its echo integration window is not stated (this example integrates the complex echo over 1.0 us after the second pulse), and its carrier step is not stated. Relaxation (T1=100 us, T2=4 us in the paper) is omitted here; with T2 four times the 1.0 us echo window it scales the four spectra by nearly the same factor, an approximation rather than an identity.

## Numerical / algorithmic content

- Time propagation is explicit in Hilbert space (`zeeman-hilb`): the pulse sequence `experiments/echo_sweep.m` steps through the Hamiltonian rotor stack that `singlerot()` supplies, applying a propagator per 5 ns time step to the density matrix.
- The rotor stack index at each time step follows from `parameters.rate`, so `rate=0` is the static case; the rotor rank (2700) is set so that the rotor phase resolution of the stack matches the 5 ns time step at 37 kHz, the fastest rate used; at slower rates consecutive steps reuse a stack element.
- The electron coherence pathway (-1 after the first pulse, +1 after the second) is selected with `coherence()` in place of the phase cycle; the sequence averages over 100 rotor phases at the start of the sequence and sweeps the carrier offset inside the sequence, reusing the free-evolution propagators because the rotor stack commutes with the electron `Lz`.

## Implementation structure

- Two-pulse echo-detected frequency-swept EPR spectra of the P1 sub-
- stitutional nitrogen defect in diamond, static and under magic ang-
- le spinning, after Figure 1a of Khamrui et al., J. Phys. Chem. Lett.
- 2026, <https://doi.org/10.1021/acs.jpclett.6c02108>: 400 ns long
- pulses with 416 kHz nutation frequency, 300 ns interpulse delay, at
- 6.9 T, static and at 10, 25, and 37 kHz MAS. The carrier is stepped
- across the spectrum and the integrated echo is recorded at each fre-
- quency, with the spectra normalised to the static one.
- The 14N hyperfine coupling (dipolar part 10.9 MHz) makes the two outer
- lines dephase under spinning as their resonance frequencies move during
- the sequence, whereas the central line survives. Compared to the paper's
- own simulation, the static outer edges here are weaker (0.2 of the cent-
- ral peak against 0.4) and the central line dephases less (0.8 at 37 kHz).
- Calculation time: hours on a 256-core node.
- P1 centre parameters
- Build the spin system
- Magnet field, central line at 193.797 GHz
- Basis set
- Spinach housekeeping
- Rotor parameters
- Sequence parameters
- Spinning rates
- Simulation
- Normalisation to the static spectrum
- Plotting
