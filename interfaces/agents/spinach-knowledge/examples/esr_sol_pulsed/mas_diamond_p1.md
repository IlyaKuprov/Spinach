# examples/esr_sol_pulsed/mas_diamond_p1.m

- Signature: `mas_diamond_p1()`

## Purpose

Calculates two-pulse echo-detected, frequency-swept EPR spectra of the P1 substitutional nitrogen defect in diamond, both static and under magic-angle spinning (MAS). It follows Figure 1a of Khamrui et al., *J. Phys. Chem. Lett.* (2026), [doi:10.1021/acs.jpclett.6c02108](https://doi.org/10.1021/acs.jpclett.6c02108). The example uses 400 ns pulses with a 416 kHz nutation frequency, a 300 ns interpulse delay, and a 6.9156 T field; the source identifies the central line at 193.797 GHz. MAS rates are 10, 25, and 37 kHz. The carrier is swept and the integrated echo is recorded at each offset. Each spectrum is normalised by the maximum of the static spectrum. Calculation time: hours on a 256-core node.

## Model and sequence

The P1 system is created for a `14N` centre with orientation `111`. The example uses the full Zeeman Hilbert-space basis. It simulates rates [0, 10, 25, 37] kHz with the rotor axis along [1, 1, 1], rotor-stack rank 2700, and 100 initial rotor phases. The sequence has a 300 ns interpulse delay and integrates the echo for 1.0 μs after the second pulse, using 5 ns time steps. A 400 ns pulse has a 416 kHz nutation frequency. The carrier sweep is 300 MHz, sampled at 601 points on a GHz lab-frame axis.

The Hamiltonian rotor stack from `singlerot` is advanced by `echo_sweep`. It selects the electron coherence pathway (−1 after the first pulse, +1 after the second) instead of using a phase cycle; the axial, coaxial P1 tensors permit a two-angle powder grid. Relaxation is omitted.

## Comparison with Figure 1a

The line positions and collapse of the outer lines under spinning match the paper. The 14N hyperfine coupling, with dipolar part 10.9 MHz, makes the outer lines dephase as their resonance frequencies move during the sequence, while the central line survives. However, the intensities do not reproduce the paper's simulation. Relative to the static central peak, the static outer perpendicular edges are 0.20 and 0.18 here versus about 0.4 in the paper. The central line retains 0.98, 0.90, and 0.82 of its echo at 10, 25, and 37 kHz, versus about 0.87, 0.52, and 0.33 in the paper. The paper's simulated outer lines are broad humps; this example gives a perpendicular-edge singularity convolved with the roughly 1 MHz pulse response.

Zeroing the 14N quadrupole leaves the 37 kHz central-line survival unchanged, so the quadrupole is not the source of that discrepancy. Other candidates remain untested: the paper's model retains only the secular hyperfine term and omits nuclear Zeeman interaction; its echo-integration window and carrier step are not stated. The paper uses T1 = 100 μs and T2 = 4 μs, while this example omits relaxation. Because T2 is four times the 1.0 μs echo window, treating relaxation as a common factor across the four spectra is an approximation, not an identity.
