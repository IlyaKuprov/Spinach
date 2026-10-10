# examples/esr_sol_pulsed/mas_diamond_p1.m

- Signature: `mas_diamond_p1()`
- Source: [`examples/esr_sol_pulsed/mas_diamond_p1.m`](https://github.com/IlyaKuprov/Spinach/blob/main/examples/esr_sol_pulsed/mas_diamond_p1.m)
- Spin-system builder: [`etc/diamond_defects/diamond_p1.m`](https://github.com/IlyaKuprov/Spinach/blob/main/etc/diamond_defects/diamond_p1.m)
- Sequence: [`experiments/echo_sweep.m`](https://github.com/IlyaKuprov/Spinach/blob/main/experiments/echo_sweep.m)
- Figure comparison: Khamrui et al., *J. Phys. Chem. Lett.* (2026), [doi:10.1021/acs.jpclett.6c02108](https://doi.org/10.1021/acs.jpclett.6c02108). The spin-system builder cites magnetic parameters from [Nir-Arad et al. (2024)](https://doi.org/10.1039/d4cp03055a) and [Smith et al. (1959)](https://doi.org/10.1103/PhysRev.115.1546).

## Physical aim and spin model

The calculation compares two-pulse, echo-detected frequency-swept EPR of a diamond P1 substitutional-nitrogen centre at rest and under magic-angle spinning. The `diamond_p1` builder supplies an electron and `14N`, with the P1 orientation set to `111`. Its electron g principal values are [2.00220, 2.00220, 2.00218]; the nitrogen hyperfine tensor principal values are [81.3, 81.3, 114.0] MHz and the nitrogen quadrupole term is built from `D = −3.97 MHz`, `E = 0`. The tensors are axial/coaxial in the builder's crystal frame. The source describes the dipolar part of the nitrogen hyperfine coupling as 10.9 MHz. The full Zeeman Hilbert-space basis is used without approximation.

## Echo sweep and plotted result

The field is 6.9156 T (the source identifies the central line at 193.797 GHz). Two 400 ns pulses, with 416 kHz nutation frequency, are separated by 300 ns; the echo is integrated for 1.0 μs after the second pulse in 5 ns time steps. The rotor axis is [1,1,1], rotor-stack maximum rank is 2700, and 100 initial rotor phases are averaged. The carrier is swept over 300 MHz in 601 points, with zero offset and the `rep_2ang_400pts_sph` grid. The `singlerot`/`echo_sweep` path computes the integrated echo at each carrier position; the electron coherence pathway (−1 then +1) replaces phase cycling. The rates are 0, 10, 25, and 37 kHz. Absolute spectra are normalised to the maximum of the static spectrum, then plotted together with a legend and echo-intensity axis in arbitrary units.

Relaxation is omitted, and the plot is the only output (no spectrum file is saved). The source header's comparison describes static outer edges at 0.2 of the central peak versus about 0.4 in the paper, and about 0.8 of the central-line echo remaining at 37 kHz; it notes that the paper's simulated result differs. The header estimates hours on a 256-core node.
